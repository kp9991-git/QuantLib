/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2021 Marcin Rybacki
 Copyright (C) 2025 Uzair Beg
 Copyright (C) 2026 Kyrylo Protsenko

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/cashflows/couponsensitivities.hpp>
#include <ql/cashflows/floatingratecoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/currencies/exchangeratemanager.hpp>
#include <ql/experimental/fx/discountingmtmcrosscurrencybasisswapengine.hpp>
#include <ql/experimental/fx/fxresetcashflows.hpp>
#include <ql/experimental/termstructures/crosscurrencyratehelpers.hpp>
#include <ql/money.hpp>
#include <ql/pricingengines/swap/discountingconstnotionalcrosscurrencyswapengine.hpp>
#include <ql/utilities/null_deleter.hpp>
#include <utility>

namespace QuantLib {

    namespace {

        constexpr double sample_fixed_rate = 0.01;

        // Treat an explicitly-passed NoFrequency the same as an unset (nullopt)
        // payment frequency.  Before these parameters were migrated to
        // std::optional<Frequency>, NoFrequency was the sentinel meaning "derive
        // the schedule from the index tenor".  Normalizing it here preserves that
        // behavior and keeps the stored optional either empty or holding an
        // actual frequency, so the rest of the code can treat the two cases
        // identically.
        std::optional<Frequency> normalizedPaymentFrequency(std::optional<Frequency> frequency) {
            if (frequency && *frequency == NoFrequency)
                return std::nullopt;
            return frequency;
        }

        Schedule legSchedule(const Date& evaluationDate,
                             const Period& tenor,
                             const Period& frequency,
                             Natural fixingDays,
                             const Calendar& calendar,
                             BusinessDayConvention convention,
                             bool endOfMonth) {
            QL_REQUIRE(tenor >= frequency,
                       "XCCY instrument tenor should not be smaller than coupon frequency.");

            Date referenceDate = calendar.adjust(evaluationDate);
            Date earliestDate = calendar.advance(referenceDate, fixingDays * Days, convention);
            Date maturity = earliestDate + tenor;
            return MakeSchedule()
                .from(earliestDate)
                .to(maturity)
                .withTenor(frequency)
                .withCalendar(calendar)
                .withConvention(convention)
                .endOfMonth(endOfMonth)
                .backwards();
        }

        Schedule floatingLegSchedule(const Date& evaluationDate,
                                     const Period& tenor,
                                     Natural fixingDays,
                                     const Calendar& calendar,
                                     BusinessDayConvention convention,
                                     bool endOfMonth,
                                     const ext::shared_ptr<IborIndex>& idx,
                                     std::optional<Frequency> paymentFrequency) {
            auto overnightIndex = ext::dynamic_pointer_cast<OvernightIndex>(idx);

            Period freqPeriod;
            if (!paymentFrequency) {
                QL_REQUIRE(!overnightIndex, "Require payment frequency for overnight indices.");
                freqPeriod = idx->tenor();
            } else {
                freqPeriod = Period(*paymentFrequency);
            }

            return legSchedule(evaluationDate, tenor, freqPeriod, fixingDays, calendar,
                               convention, endOfMonth);
        }

        Leg buildFloatingLeg(const Schedule& schedule,
                             const ext::shared_ptr<IborIndex>& idx,
                             Integer paymentLag,
                             std::optional<bool> useIndexedCoupons) {
            if (auto overnightIndex = ext::dynamic_pointer_cast<OvernightIndex>(idx)) {
                return OvernightLeg(schedule, overnightIndex)
                    .withNotionals(1.0)
                    .withPaymentLag(paymentLag);
            }
            return IborLeg(schedule, idx)
                .withNotionals(1.0)
                .withPaymentLag(paymentLag)
                .withIndexedCoupons(useIndexedCoupons);
        }

        std::pair<Real, Real>
        npvbpsConstNotionalLeg(const Leg& leg,
                               const Date& initialNotionalExchangeDate,
                               const Date& finalNotionalExchangeDate,
                               const Handle<YieldTermStructure>& discountCurveHandle) {
            const Spread basisPoint = 1.0e-4;
            Date refDt = discountCurveHandle->referenceDate();
            const YieldTermStructure& discountRef = **discountCurveHandle;
            bool includeSettleDtFlows = true;
            auto [npv, bps] = CashFlows::npvbps(leg, discountRef, includeSettleDtFlows, refDt, refDt);
            // Include NPV of the notional exchange at start and maturity.
            npv += (-1.0) * discountRef.discount(initialNotionalExchangeDate);
            npv += discountRef.discount(finalNotionalExchangeDate);
            bps /= basisPoint;
            return { npv, bps };
        }

    }


    CrossCurrencySwapRateHelperBase::CrossCurrencySwapRateHelperBase(
        const Handle<Quote>& quote,
        const Period& tenor,
        Natural fixingDays,
        Calendar calendar,
        BusinessDayConvention convention,
        bool endOfMonth,
        Handle<YieldTermStructure> collateralCurve,
        Integer paymentLag,
        bool paymentLagOnNotionalExchanges)
    : RelativeDateRateHelper(quote), tenor_(tenor), fixingDays_(fixingDays),
      calendar_(std::move(calendar)), convention_(convention), endOfMonth_(endOfMonth),
      paymentLag_(paymentLag), paymentLagOnNotionalExchanges_(paymentLagOnNotionalExchanges),
      collateralHandle_(std::move(collateralCurve)) {
        registerWith(collateralHandle_);
    }

    void CrossCurrencySwapRateHelperBase::setTermStructure(YieldTermStructure* t) {
        // do not set the relinkable handle as an observer -
        // force recalculation when needed
        bool observer = false;

        ext::shared_ptr<YieldTermStructure> temp(t, null_deleter());
        termStructureHandle_.linkTo(temp, observer);

        RelativeDateRateHelper::setTermStructure(t);
    }

    void CrossCurrencySwapRateHelperBase::initializeDatesFromLegs(const Leg& firstLeg,
                                                                  const Leg& secondLeg) {
        earliestDate_ = std::min(CashFlows::startDate(firstLeg),
                                 CashFlows::startDate(secondLeg));

        maturityDate_ = std::max(CashFlows::maturityDate(firstLeg),
                                 CashFlows::maturityDate(secondLeg));

        // Principal exchanges settle on the effective and maturity dates:
        // the payment lag applies only to coupons -- unless the swap's
        // convention lags the principal flows too, in which case both
        // exchanges move by the same lag as the coupons and the final
        // exchange settles together with the final coupon.  advance() with
        // zero days reduces to adjust(), preserving the default behaviour.
        Integer exchangeLag = paymentLagOnNotionalExchanges_ ? paymentLag_ : 0;
        initialNotionalExchangeDate_ = calendar_.advance(earliestDate_, exchangeLag, Days, convention_);
        finalNotionalExchangeDate_   = calendar_.advance(maturityDate_, exchangeLag, Days, convention_);

        Date lastPaymentDate =
            std::max(firstLeg.back()->date(),
                     secondLeg.back()->date());

        latestRelevantDate_ = latestDate_ = std::max(maturityDate_, lastPaymentDate);
    }



    CrossCurrencyBasisSwapRateHelperBase::CrossCurrencyBasisSwapRateHelperBase(
        const Handle<Quote>& basis,
        const Period& tenor,
        Natural fixingDays,
        Calendar calendar,
        BusinessDayConvention convention,
        bool endOfMonth,
        ext::shared_ptr<IborIndex> baseCurrencyIndex,
        ext::shared_ptr<IborIndex> quoteCurrencyIndex,
        Handle<YieldTermStructure> collateralCurve,
        bool isFxBaseCurrencyCollateralCurrency,
        bool isBasisOnFxBaseCurrencyLeg,
        std::optional<Frequency> paymentFrequency,
        Integer paymentLag,
        std::optional<Frequency> quoteCurrencyPaymentFrequency,
        std::optional<bool> useIndexedCoupons,
        bool paymentLagOnNotionalExchanges)
    : CrossCurrencySwapRateHelperBase(basis, tenor, fixingDays, std::move(calendar), convention, endOfMonth,
                                      std::move(collateralCurve), paymentLag,
                                      paymentLagOnNotionalExchanges),
      baseCcyIdx_(std::move(baseCurrencyIndex)), quoteCcyIdx_(std::move(quoteCurrencyIndex)),
      isFxBaseCurrencyCollateralCurrency_(isFxBaseCurrencyCollateralCurrency),
      isBasisOnFxBaseCurrencyLeg_(isBasisOnFxBaseCurrencyLeg),
      paymentFrequency_(normalizedPaymentFrequency(paymentFrequency)),
      quoteCcyPaymentFrequency_(normalizedPaymentFrequency(quoteCurrencyPaymentFrequency)),
      useIndexedCoupons_(useIndexedCoupons) {
        registerWith(baseCcyIdx_);
        registerWith(quoteCcyIdx_);

        CrossCurrencyBasisSwapRateHelperBase::initializeDates();
    }

    void CrossCurrencyBasisSwapRateHelperBase::initializeDates() {
        baseCcySchedule_ = floatingLegSchedule(evaluationDate_, tenor_, fixingDays_, calendar_,
                                               convention_, endOfMonth_, baseCcyIdx_,
                                               paymentFrequency_);
        baseCcyIborLeg_ = buildFloatingLeg(
            baseCcySchedule_, baseCcyIdx_, paymentLag_, useIndexedCoupons_);

        // If no quote-currency payment frequency was given, fall back to the
        // base-currency payment frequency (which may itself be unset, in which
        // case the quote-currency leg uses its own index tenor).
        std::optional<Frequency> effectiveQuoteCcyFreq =
            quoteCcyPaymentFrequency_ ? quoteCcyPaymentFrequency_ : paymentFrequency_;
        quoteCcySchedule_ = floatingLegSchedule(evaluationDate_, tenor_, fixingDays_, calendar_,
                                                convention_, endOfMonth_, quoteCcyIdx_,
                                                effectiveQuoteCcyFreq);
        quoteCcyIborLeg_ = buildFloatingLeg(
            quoteCcySchedule_, quoteCcyIdx_, paymentLag_, useIndexedCoupons_);

        initializeDatesFromLegs(baseCcyIborLeg_, quoteCcyIborLeg_);
    }

    const Handle<YieldTermStructure>&
    CrossCurrencyBasisSwapRateHelperBase::baseCcyLegDiscountHandle() const {
        return isFxBaseCurrencyCollateralCurrency_ ? collateralHandle_ : termStructureHandle_;
    }

    const Handle<YieldTermStructure>&
    CrossCurrencyBasisSwapRateHelperBase::quoteCcyLegDiscountHandle() const {
        return isFxBaseCurrencyCollateralCurrency_ ? termStructureHandle_ : collateralHandle_;
    }

    ConstNotionalCrossCurrencyBasisSwapRateHelper::ConstNotionalCrossCurrencyBasisSwapRateHelper(
        const Handle<Quote>& basis,
        const Period& tenor,
        Natural fixingDays,
        const Calendar& calendar,
        BusinessDayConvention convention,
        bool endOfMonth,
        const ext::shared_ptr<IborIndex>& baseCurrencyIndex,
        const ext::shared_ptr<IborIndex>& quoteCurrencyIndex,
        const Handle<YieldTermStructure>& collateralCurve,
        bool isFxBaseCurrencyCollateralCurrency,
        bool isBasisOnFxBaseCurrencyLeg,
        std::optional<Frequency> paymentFrequency,
        Integer paymentLag,
        std::optional<Frequency> quoteCurrencyPaymentFrequency,
        std::optional<bool> useIndexedCoupons,
        bool paymentLagOnNotionalExchanges)
    : CrossCurrencyBasisSwapRateHelperBase(basis,
                                           tenor,
                                           fixingDays,
                                           calendar,
                                           convention,
                                           endOfMonth,
                                           baseCurrencyIndex,
                                           quoteCurrencyIndex,
                                           collateralCurve,
                                           isFxBaseCurrencyCollateralCurrency,
                                           isBasisOnFxBaseCurrencyLeg,
                                           paymentFrequency,
                                           paymentLag,
                                           quoteCurrencyPaymentFrequency,
                                           useIndexedCoupons,
                                           paymentLagOnNotionalExchanges) {
        buildSwap();
    }

    void ConstNotionalCrossCurrencyBasisSwapRateHelper::initializeDates() {
        CrossCurrencyBasisSwapRateHelperBase::initializeDates();
        buildSwap();
    }

    void ConstNotionalCrossCurrencyBasisSwapRateHelper::buildSwap() {
        // The exposed swap mirrors the helper's par convention: unit notionals,
        // zero spreads and spot FX = 1, so that its fair spread on the basis
        // leg reproduces the helper quote.  It pays the base-currency leg.
        swap_ = ext::make_shared<ConstNotionalCrossCurrencyBasisSwap>(
            1.0, baseCcyIdx_->currency(), baseCcySchedule_, baseCcyIdx_, 0.0, 1.0,
            1.0, quoteCcyIdx_->currency(), quoteCcySchedule_, quoteCcyIdx_, 0.0, 1.0,
            paymentLag_, paymentLag_, false, Null<Natural>(), false, 0,
            RateAveraging::Compound, false, Null<Natural>(), false, 0,
            RateAveraging::Compound, false, useIndexedCoupons_,
            paymentLagOnNotionalExchanges_);
        swap_->setPricingEngine(ext::make_shared<DiscountingConstNotionalCrossCurrencySwapEngine>(
            quoteCcyIdx_->currency(), quoteCcyLegDiscountHandle(),
            baseCcyIdx_->currency(), baseCcyLegDiscountHandle(),
            makeQuoteHandle(1.0), true));
    }

    Real ConstNotionalCrossCurrencyBasisSwapRateHelper::impliedQuote() const {
        QL_REQUIRE(!termStructureHandle_.empty(), "term structure not set");
        QL_REQUIRE(!collateralHandle_.empty(), "collateral term structure not set");

        auto [npvBaseCcy, bpsBaseCcy] = npvbpsConstNotionalLeg(baseCcyIborLeg_, initialNotionalExchangeDate_, finalNotionalExchangeDate_, baseCcyLegDiscountHandle());
        auto [npvQuoteCcy, bpsQuoteCcy] = npvbpsConstNotionalLeg(quoteCcyIborLeg_, initialNotionalExchangeDate_, finalNotionalExchangeDate_, quoteCcyLegDiscountHandle());

        Real bps = isBasisOnFxBaseCurrencyLeg_ ? -bpsBaseCcy : bpsQuoteCcy;

        QL_REQUIRE(std::fabs(bps) > 0.0, "null BPS");

        return -(npvQuoteCcy - npvBaseCcy) / bps;
    }

    QuoteSensitivities
    ConstNotionalCrossCurrencyBasisSwapRateHelper::impliedQuoteSensitivitiesByCurve() const {
        if (termStructure_ == nullptr || termStructureHandle_.empty() ||
            collateralHandle_.empty())
            return {};

        // Q = -(npvQuote - npvBase)/bps
        // bps is -A_base or A_quote, depending on the quoted leg
        QuoteSensitivities result;
        const Leg* legs[2] = {&baseCcyIborLeg_, &quoteCcyIborLeg_};
        const Handle<YieldTermStructure> discountHandles[2] = {
            baseCcyLegDiscountHandle(), quoteCcyLegDiscountHandle()};

        std::vector<detail::FloatingFlowData> flows[2];
        for (Size legNo = 0; legNo < 2; ++legNo) {
            const YieldTermStructure& discountCurve = **discountHandles[legNo];
            if (!detail::analyzeFloatingLeg(*legs[legNo], discountCurve.referenceDate(),
                                            discountCurve, flows[legNo], result,
                                            /*includeSettlementDateFlows=*/true))
                return {};
        }

        Size basisLegNo = isBasisOnFxBaseCurrencyLeg_ ? 0 : 1;
        Real annuitySign = basisLegNo == 0 ? -1.0 : 1.0;
        Real annuity = 0.0;
        for (const auto& d : flows[basisLegNo])
            annuity += d.ntau*d.discount;
        Real bps = annuitySign*annuity;
        if (bps == 0.0)
            return {};

        Real quote = impliedQuote();
        // dQ/dP = -d(npvQuote - npvBase)/dP / bps - Q * dbps/dP / bps
        for (Size legNo = 0; legNo < 2; ++legNo) {
            Real npvSign = legNo == 0 ? 1.0 : -1.0;
            auto& discountBucket = result.sensitivities[&**discountHandles[legNo]];
            for (const auto& d : flows[legNo]) {
                Real derivative = npvSign*d.amount/bps;
                if (legNo == basisLegNo)
                    derivative -= quote*annuitySign*d.ntau/bps;
                discountBucket.emplace_back(d.payDate, derivative);
            }
            // notional exchanges at start and maturity
            discountBucket.emplace_back(initialNotionalExchangeDate_, -npvSign/bps);
            discountBucket.emplace_back(finalNotionalExchangeDate_, npvSign/bps);
            // coupon forecast sensitivities
            for (const auto& d : flows[legNo]) {
                if (d.sensitivities.empty())
                    continue;
                auto& bucket = result.sensitivities[
                    static_cast<const TermStructure*>(d.forecastCurve)];
                for (const auto& [date, w] : d.sensitivities)
                    bucket.emplace_back(date, npvSign*d.discount*w/bps);
            }
        }
        result.available = true;
        return result;
    }

    void ConstNotionalCrossCurrencyBasisSwapRateHelper::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<ConstNotionalCrossCurrencyBasisSwapRateHelper>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            RateHelper::accept(v);
    }


    MtMCrossCurrencyBasisSwapRateHelper::MtMCrossCurrencyBasisSwapRateHelper(
        const Handle<Quote>& basis,
        const Period& tenor,
        Natural fixingDays,
        const Calendar& calendar,
        BusinessDayConvention convention,
        bool endOfMonth,
        const ext::shared_ptr<IborIndex>& baseCurrencyIndex,
        const ext::shared_ptr<IborIndex>& quoteCurrencyIndex,
        const Handle<YieldTermStructure>& collateralCurve,
        bool isFxBaseCurrencyCollateralCurrency,
        bool isBasisOnFxBaseCurrencyLeg,
        bool isFxBaseCurrencyLegResettable,
        std::optional<Frequency> paymentFrequency,
        Integer paymentLag,
        std::optional<Frequency> quoteCurrencyPaymentFrequency,
        Natural fxResetFixingDays,
        Calendar fxResetFixingCalendar,
        std::optional<bool> useIndexedCoupons)
    : CrossCurrencyBasisSwapRateHelperBase(basis,
                                           tenor,
                                           fixingDays,
                                           calendar,
                                           convention,
                                           endOfMonth,
                                           baseCurrencyIndex,
                                           quoteCurrencyIndex,
                                           collateralCurve,
                                           isFxBaseCurrencyCollateralCurrency,
                                           isBasisOnFxBaseCurrencyLeg,
                                           paymentFrequency,
                                           paymentLag,
                                           quoteCurrencyPaymentFrequency,
                                           useIndexedCoupons),
      isFxBaseCurrencyLegResettable_(isFxBaseCurrencyLegResettable),
      fxResetFixingDays_(fxResetFixingDays), fxResetFixingCalendar_(std::move(fxResetFixingCalendar)) {
        buildSwap();
    }

    void MtMCrossCurrencyBasisSwapRateHelper::initializeDates() {
        CrossCurrencyBasisSwapRateHelperBase::initializeDates();
        buildSwap();
    }

    void MtMCrossCurrencyBasisSwapRateHelper::buildSwap() {
        // The exposed swap mirrors the helper's par convention: unit notionals,
        // zero spreads and spot FX = 1, so that its fair spread on the basis
        // leg reproduces the helper quote.  It pays the base-currency leg.
        swap_ = ext::make_shared<MtMCrossCurrencyBasisSwap>(
            MtMCrossCurrencyBasisSwap::Type::PayFxBaseCurrency,
            1.0, baseCcyIdx_->currency(), baseCcySchedule_, baseCcyIdx_, 0.0, 1.0,
            1.0, quoteCcyIdx_->currency(), quoteCcySchedule_, quoteCcyIdx_, 0.0, 1.0,
            isFxBaseCurrencyLegResettable_, fxResetFixingDays_,
            fxResetFixingCalendar_, paymentLag_, paymentLag_,
            convention_, convention_, false, Null<Natural>(), false, 0,
            RateAveraging::Compound, false, Null<Natural>(), false, 0,
            RateAveraging::Compound, false, useIndexedCoupons_);
        swap_->setPricingEngine(ext::make_shared<DiscountingMtMCrossCurrencyBasisSwapEngine>(
            quoteCcyIdx_->currency(), quoteCcyLegDiscountHandle(),
            baseCcyIdx_->currency(), baseCcyLegDiscountHandle(),
            makeQuoteHandle(1.0), true));
    }

    Real MtMCrossCurrencyBasisSwapRateHelper::impliedQuote() const {
        QL_REQUIRE(!termStructureHandle_.empty(), "term structure not set");
        QL_REQUIRE(!collateralHandle_.empty(), "collateral term structure not set");

        swap_->deepUpdate();
        if (isBasisOnFxBaseCurrencyLeg_)
            return swap_->fairFxBaseSpread();
        return swap_->fairFxQuoteSpread();
    }

    QuoteSensitivities
    MtMCrossCurrencyBasisSwapRateHelper::impliedQuoteSensitivitiesByCurve() const {
        if (termStructure_ == nullptr || termStructureHandle_.empty() ||
            collateralHandle_.empty())
            return {};

        // Both discount curves drive NPVs, FX conversion, and reset notionals
        // Coupon forwards add the index forecast curves
        // Q = -NPV/B with
        //   NPV = -fx*N_base + N_quote,   B = payer_k*fxconv_k*A_k
        QuoteSensitivities result;
        const Handle<YieldTermStructure>& baseDisc = baseCcyLegDiscountHandle();
        const Handle<YieldTermStructure>& quoteDisc = quoteCcyLegDiscountHandle();
        const auto* baseKey = static_cast<const TermStructure*>(&**baseDisc);
        const auto* quoteKey = static_cast<const TermStructure*>(&**quoteDisc);
        const YieldTermStructure& curve = **termStructureHandle_;
        Date refDate = curve.referenceDate();
        Date today = Settings::instance().evaluationDate();

        Size resettingLegNo = isFxBaseCurrencyLegResettable_ ? 0 : 1;
        const Handle<YieldTermStructure>& resetCurve =
            resettingLegNo == 0 ? baseDisc : quoteDisc;
        const Handle<YieldTermStructure>& constCurve =
            resettingLegNo == 0 ? quoteDisc : baseDisc;
        const auto* resetKey = static_cast<const TermStructure*>(&**resetCurve);
        const auto* constKey = static_cast<const TermStructure*>(&**constCurve);
        const Currency& resetCcy =
            resettingLegNo == 0 ? baseCcyIdx_->currency() : quoteCcyIdx_->currency();
        const Currency& constCcy =
            resettingLegNo == 0 ? quoteCcyIdx_->currency() : baseCcyIdx_->currency();

        // FX settlement date used by the swap engine
        Calendar fxCalendar = (fxResetFixingDays_ != 0 && fxResetFixingCalendar_.empty())
            ? calendar_ : fxResetFixingCalendar_;
        Date fxSettle = fxCalendar.empty() ? refDate :
            FxResetConvention(fxResetFixingDays_, fxCalendar).valueDate(today);

        using TaggedEntry = std::tuple<const TermStructure*, Date, Real>;

        // base-to-quote conversion at unit spot
        Real fx = quoteDisc->discount(fxSettle)/baseDisc->discount(fxSettle);
        const TaggedEntry dFx[2] = {
            {quoteKey, fxSettle, fx/quoteDisc->discount(fxSettle)},
            {baseKey, fxSettle, -fx/baseDisc->discount(fxSettle)}};

        // same convention as DiscountingFxResetPricer::fxRate at unit spot
        auto historicalReset = [&](const Date& fixingDate) -> Real {
            try {
                ExchangeRate rate = ExchangeRateManager::instance().lookup(
                    constCcy, resetCcy, fixingDate);
                return rate.exchange(Money(1.0, constCcy)).value();
            } catch (Error&) {
                return Null<Real>();
            }
        };
        struct Reset { Real rate; bool forecast; };
        auto analyzeReset = [&](const FxReset& reset) -> Reset {
            if (reset.fixingDate() <= today) {
                Real r = historicalReset(reset.fixingDate());
                if (r != Null<Real>())
                    return {r, false};
                if (reset.fixingDate() < today ||
                    Settings::instance().enforcesTodaysHistoricFixings())
                    return {Null<Real>(), false};  // required fixing is missing
            }
            Real rate = (resetCurve->discount(fxSettle)/constCurve->discount(fxSettle))
                * (constCurve->discount(reset.valueDate())
                   /resetCurve->discount(reset.valueDate()));
            return {rate, true};
        };
        // d(scale * fxRate)/dP for both discount curves
        auto addResetSensitivities = [&](const FxReset& reset, const Reset& r, Real scale,
                                         std::vector<TaggedEntry>& out) {
            if (!r.forecast)
                return;
            out.emplace_back(resetKey, fxSettle,
                             scale*r.rate/resetCurve->discount(fxSettle));
            out.emplace_back(resetKey, reset.valueDate(),
                             -scale*r.rate/resetCurve->discount(reset.valueDate()));
            out.emplace_back(constKey, fxSettle,
                             -scale*r.rate/constCurve->discount(fxSettle));
            out.emplace_back(constKey, reset.valueDate(),
                             scale*r.rate/constCurve->discount(reset.valueDate()));
        };

        struct FlowData {
            Date payDate;
            Real amount = 0.0, ntau = 0.0;
            std::vector<TaggedEntry> dAmount, dNtau;
        };
        // fall back to amount() and mark incomplete forecast curves
        auto analyzeFlowCoupon = [&](const ext::shared_ptr<CashFlow>& cf,
                                     detail::CouponSensitivityAnalysis& a) -> bool {
            a = detail::analyzeCoupon(cf, true);
            if (!a.supported) {
                a = detail::analyzeCoupon(cf, false);
                if (!a.supported)
                    return false;
                if (ext::dynamic_pointer_cast<FloatingRateCoupon>(cf) != nullptr) {
                    if (a.forecastCurve == nullptr)
                        // unidentified forecast curve
                        return false;
                    result.incomplete.insert(a.forecastCurve);
                }
                a.amountSensitivities.clear();
            }
            return true;
        };

        Real legNPV[2] = {0.0, 0.0}, legAnnuity[2] = {0.0, 0.0};
        std::vector<FlowData> flows[2];
        for (Size legNo = 0; legNo < 2; ++legNo) {
            const Handle<YieldTermStructure>& legCurve = legNo == 0 ? baseDisc : quoteDisc;
            for (const auto& cf : swap_->leg(legNo)) {
                if (cf->hasOccurred(refDate, true))
                    continue;
                FlowData d;
                d.payDate = cf->date();
                if (auto coupon = ext::dynamic_pointer_cast<FxResetCoupon>(cf)) {
                    auto r = analyzeReset(coupon->fxReset());
                    if (r.rate == Null<Real>())
                        return {};
                    // rescale the underlying coupon to the reset notional
                    detail::CouponSensitivityAnalysis ua;
                    if (!analyzeFlowCoupon(coupon->underlying(), ua))
                        return {};
                    Real underlyingScale =
                        coupon->constantLegNotional()/coupon->underlying()->nominal();
                    Real underlyingAmount = ua.amount;
                    d.amount = underlyingAmount*underlyingScale*r.rate;
                    d.ntau = coupon->constantLegNotional()*r.rate*coupon->accrualPeriod();
                    addResetSensitivities(coupon->fxReset(), r,
                                          underlyingAmount*underlyingScale, d.dAmount);
                    addResetSensitivities(coupon->fxReset(), r,
                                          coupon->constantLegNotional()*coupon->accrualPeriod(),
                                          d.dNtau);
                    // underlying coupon forecast sensitivity
                    for (const auto& [date, w] : ua.amountSensitivities)
                        d.dAmount.emplace_back(
                            static_cast<const TermStructure*>(ua.forecastCurve),
                            date, underlyingScale*r.rate*w);
                } else if (auto exchange =
                               ext::dynamic_pointer_cast<FxResetNotionalExchange>(cf)) {
                    if (exchange->previousReset()) {
                        auto r = analyzeReset(*exchange->previousReset());
                        if (r.rate == Null<Real>())
                            return {};
                        d.amount += exchange->constantLegNotional()*r.rate;
                        addResetSensitivities(*exchange->previousReset(), r,
                                              exchange->constantLegNotional(), d.dAmount);
                    }
                    if (exchange->currentReset()) {
                        auto r = analyzeReset(*exchange->currentReset());
                        if (r.rate == Null<Real>())
                            return {};
                        d.amount -= exchange->constantLegNotional()*r.rate;
                        addResetSensitivities(*exchange->currentReset(), r,
                                              -exchange->constantLegNotional(), d.dAmount);
                    }
                } else if (ext::dynamic_pointer_cast<Coupon>(cf) != nullptr) {
                    detail::CouponSensitivityAnalysis a;
                    if (!analyzeFlowCoupon(cf, a))
                        return {};
                    d.amount = a.amount;
                    d.ntau = a.ntau;
                    for (const auto& [date, w] : a.amountSensitivities)
                        d.dAmount.emplace_back(
                            static_cast<const TermStructure*>(a.forecastCurve), date, w);
                } else {
                    d.amount = cf->amount();
                }
                DiscountFactor P = legCurve->discount(d.payDate);
                legNPV[legNo] += d.amount*P;
                legAnnuity[legNo] += d.ntau*P;
                flows[legNo].push_back(std::move(d));
            }
        }

        Size basisLegNo = isBasisOnFxBaseCurrencyLeg_ ? 0 : 1;
        Real payerBasis = basisLegNo == 0 ? -1.0 : 1.0;
        Real fxconvBasis = basisLegNo == 0 ? fx : 1.0;
        Real B = payerBasis*fxconvBasis*legAnnuity[basisLegNo];
        if (B == 0.0)
            return {};
        Real quote = impliedQuote();

        // dQ = -dNPV/B - Q*dB/B
        auto emit = [&](const TermStructure* key, const Date& d, Real v) {
            result.sensitivities[key].emplace_back(d, v);
        };

        for (Size legNo = 0; legNo < 2; ++legNo) {
            const Handle<YieldTermStructure>& legCurve = legNo == 0 ? baseDisc : quoteDisc;
            const auto* legKey = legNo == 0 ? baseKey : quoteKey;
            Real payer = legNo == 0 ? -1.0 : 1.0;
            Real fxconv = legNo == 0 ? fx : 1.0;
            for (const auto& d : flows[legNo]) {
                DiscountFactor P = legCurve->discount(d.payDate);
                // amount contribution to -dNPV/B
                for (const auto& [key, date, w] : d.dAmount)
                    emit(key, date, -payer*fxconv*P*w/B);
                emit(legKey, d.payDate, -payer*fxconv*d.amount/B);
                // annuity contribution to -Q*dB/B
                if (legNo == basisLegNo) {
                    for (const auto& [key, date, w] : d.dNtau)
                        emit(key, date, -quote*payerBasis*fxconvBasis*P*w/B);
                    emit(legKey, d.payDate, -quote*payerBasis*fxconvBasis*d.ntau/B);
                }
            }
        }
        // FX-conversion contribution to NPV and base-leg B
        for (const auto& [key, date, w] : dFx)
            emit(key, date, legNPV[0]*w/B);
        if (basisLegNo == 0)
            for (const auto& [key, date, w] : dFx)
                emit(key, date, quote*legAnnuity[0]*w/B);
        result.available = true;
        return result;
    }

    void MtMCrossCurrencyBasisSwapRateHelper::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<MtMCrossCurrencyBasisSwapRateHelper>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            RateHelper::accept(v);
    }


    ConstNotionalCrossCurrencySwapRateHelper::ConstNotionalCrossCurrencySwapRateHelper(
        const Handle<Quote>& fixedRate,
        const Period& tenor,
        Natural fixingDays,
        const Calendar& calendar,
        BusinessDayConvention convention,
        bool endOfMonth,
        Frequency fixedFrequency,
        DayCounter  fixedDayCount,
        const ext::shared_ptr<IborIndex>& floatIndex,
        const Handle<YieldTermStructure>& collateralCurve,
        bool collateralOnFixedLeg,
        Integer paymentLag,
        std::optional<bool> useIndexedCoupons,
        std::optional<Frequency> floatPaymentFrequency)
    : CrossCurrencySwapRateHelperBase(fixedRate, tenor, fixingDays, calendar, convention, endOfMonth,
                                      collateralCurve, paymentLag),
      fixedFrequency_(fixedFrequency),
      fixedDayCount_(std::move(fixedDayCount)),
      floatIndex_(floatIndex),
      collateralOnFixedLeg_(collateralOnFixedLeg),
      useIndexedCoupons_(useIndexedCoupons),
      floatPaymentFrequency_(normalizedPaymentFrequency(floatPaymentFrequency)) {

        QL_REQUIRE(floatIndex_, "floating index required");
        registerWith(floatIndex_);

        initializeDates();
    }

    void ConstNotionalCrossCurrencySwapRateHelper::initializeDates() {
        Real nominal = 1.0;
        // Both legs roll on the helper's calendar and convention.  The floating
        // index supplies its fixing calendar for fixings, not for the accrual schedule.
        // Taking the roll convention or calendar from the index would make the
        // two legs diverge at a month-end roll and move the implied pillar.
        Schedule fixedSch = legSchedule(evaluationDate_, tenor_, Period(fixedFrequency_), fixingDays_, calendar_,
                                       convention_, endOfMonth_);
        // An overnight index has no payment frequency of its own, so
        // floatingLegSchedule requires one to have been supplied.
        Schedule floatSch = floatingLegSchedule(evaluationDate_, tenor_, fixingDays_, calendar_,
                                                convention_, endOfMonth_, floatIndex_,
                                                floatPaymentFrequency_);

        xccySwap_ = ext::make_shared<ConstNotionalCrossCurrencyFixedVsFloatingSwap>(
            Swap::Payer,
            nominal,
            Currency(), 
            fixedSch,
            sample_fixed_rate,
            fixedDayCount_,
            convention_,
            paymentLag_,
            calendar_,
            nominal,
            floatIndex_->currency(),
            floatSch,
            floatIndex_,
            Spread(0.0),
            convention_,
            paymentLag_,
            calendar_,
            false, false, Null<Natural>(), false, 0,
            RateAveraging::Compound, useIndexedCoupons_
        );
        auto engine = ext::make_shared<DiscountingConstNotionalCrossCurrencySwapEngine>(
            floatIndex_->currency(), floatingLegDiscountHandle(),
            Currency(), fixedLegDiscountHandle(),
            makeQuoteHandle(1.0), true);
        xccySwap_->setPricingEngine(engine);

        initializeDatesFromLegs(xccySwap_->leg(0), xccySwap_->leg(1));
    }

    const Handle<YieldTermStructure>&
    ConstNotionalCrossCurrencySwapRateHelper::fixedLegDiscountHandle() const {
        return collateralOnFixedLeg_ ? collateralHandle_ : termStructureHandle_;
    }

    const Handle<YieldTermStructure>&
    ConstNotionalCrossCurrencySwapRateHelper::floatingLegDiscountHandle() const {
        return collateralOnFixedLeg_ ? termStructureHandle_ : collateralHandle_;
    }

    Real ConstNotionalCrossCurrencySwapRateHelper::impliedQuote() const {
        QL_REQUIRE(!termStructureHandle_.empty(), "term structure not set");
        QL_REQUIRE(!collateralHandle_.empty(), "collateral term structure not set");
        xccySwap_->deepUpdate();

        return xccySwap_->fairRate();
    }

    void ConstNotionalCrossCurrencySwapRateHelper::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<ConstNotionalCrossCurrencySwapRateHelper>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            RateHelper::accept(v);
    }

}
