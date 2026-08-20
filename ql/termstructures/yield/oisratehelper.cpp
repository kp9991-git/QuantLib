/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009, 2012 Roland Lichters
 Copyright (C) 2009, 2012 Ferdinando Ametrano

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

#include <ql/instruments/makeois.hpp>
#include <ql/instruments/simplifynotificationgraph.hpp>
#include <ql/cashflows/couponpricer.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/termstructures/yield/oisratehelper.hpp>
#include <ql/utilities/null_deleter.hpp>
#include <utility>

namespace QuantLib {

    OISRateHelper::OISRateHelper(
        Natural settlementDays,
        const Period& tenor, // swap maturity
        const std::variant<Rate, Handle<Quote>>& fixedRate,
        const ext::shared_ptr<OvernightIndex>& overnightIndex,
        Handle<YieldTermStructure> discount,
        bool telescopicValueDates,
        Integer paymentLag,
        BusinessDayConvention paymentConvention,
        Frequency paymentFrequency,
        Calendar paymentCalendar,
        const Period& forwardStart,
        const std::variant<Spread, Handle<Quote>>& overnightSpread,
        Pillar::Choice pillar,
        Date customPillarDate,
        RateAveraging::Type averagingMethod,
        std::optional<bool> endOfMonth,
        std::optional<Frequency> fixedPaymentFrequency,
        Calendar fixedCalendar,
        Natural lookbackDays,
        Natural lockoutDays,
        bool applyObservationShift,
        ext::shared_ptr<FloatingRateCouponPricer> pricer,
        DateGeneration::Rule rule,
        Calendar overnightCalendar,
        BusinessDayConvention convention)
    : RelativeDateRateHelper(fixedRate), settlementDays_(settlementDays), tenor_(tenor),
      discountHandle_(std::move(discount)), telescopicValueDates_(telescopicValueDates),
      paymentLag_(paymentLag), paymentConvention_(paymentConvention),
      paymentFrequency_(paymentFrequency), paymentCalendar_(std::move(paymentCalendar)),
      forwardStart_(forwardStart), overnightSpread_(handleFromVariant(overnightSpread)), pillarChoice_(pillar),
      averagingMethod_(averagingMethod), endOfMonth_(endOfMonth),
      fixedPaymentFrequency_(fixedPaymentFrequency), fixedCalendar_(std::move(fixedCalendar)),
      overnightCalendar_(std::move(overnightCalendar)), convention_(convention),
      lookbackDays_(lookbackDays), lockoutDays_(lockoutDays),
      applyObservationShift_(applyObservationShift), pricer_(std::move(pricer)), rule_(rule) {
        initialize(overnightIndex, customPillarDate);
    }

    OISRateHelper::OISRateHelper(
        const Date& startDate,
        const Date& endDate,
        const std::variant<Rate, Handle<Quote>>& fixedRate,
        const ext::shared_ptr<OvernightIndex>& overnightIndex,
        Handle<YieldTermStructure> discount,
        bool telescopicValueDates,
        Integer paymentLag,
        BusinessDayConvention paymentConvention,
        Frequency paymentFrequency,
        Calendar paymentCalendar,
        const std::variant<Spread, Handle<Quote>>& overnightSpread,
        Pillar::Choice pillar,
        Date customPillarDate,
        RateAveraging::Type averagingMethod,
        std::optional<bool> endOfMonth,
        std::optional<Frequency> fixedPaymentFrequency,
        Calendar fixedCalendar,
        Natural lookbackDays,
        Natural lockoutDays,
        bool applyObservationShift,
        ext::shared_ptr<FloatingRateCouponPricer> pricer,
        DateGeneration::Rule rule,
        Calendar overnightCalendar,
        BusinessDayConvention convention)
    : RelativeDateRateHelper(fixedRate, false), startDate_(startDate), endDate_(endDate),
      discountHandle_(std::move(discount)), telescopicValueDates_(telescopicValueDates),
      paymentLag_(paymentLag), paymentConvention_(paymentConvention),
      paymentFrequency_(paymentFrequency), paymentCalendar_(std::move(paymentCalendar)),
      overnightSpread_(handleFromVariant(overnightSpread)), pillarChoice_(pillar),
      averagingMethod_(averagingMethod), endOfMonth_(endOfMonth),
      fixedPaymentFrequency_(fixedPaymentFrequency), fixedCalendar_(std::move(fixedCalendar)),
      overnightCalendar_(std::move(overnightCalendar)), convention_(convention),
      lookbackDays_(lookbackDays), lockoutDays_(lockoutDays),
      applyObservationShift_(applyObservationShift), pricer_(std::move(pricer)), rule_(rule) {
        initialize(overnightIndex, customPillarDate);
    }

    void OISRateHelper::initialize(const ext::shared_ptr<OvernightIndex>& overnightIndex,
                                   Date customPillarDate) {
        overnightIndex_ =
            ext::dynamic_pointer_cast<OvernightIndex>(overnightIndex->clone(termStructureHandle_));
        // We want to be notified of changes of fixings, but we don't
        // want notifications from termStructureHandle_ (they would
        // interfere with bootstrapping.)
        overnightIndex_->unregisterWith(termStructureHandle_);

        registerWith(overnightIndex_);
        registerWith(discountHandle_);
        registerWith(overnightSpread_);

        pillarDate_ = customPillarDate;
        OISRateHelper::initializeDates();
    }

    void OISRateHelper::initializeDates() {

        // 1. do not pass the spread here, as it might be a Quote
        //    i.e. it can dynamically change
        // 2. input discount curve Handle might be empty now but it could
        //    be assigned a curve later; use a RelinkableHandle here
        auto tmp = MakeOIS(tenor_, overnightIndex_, 0.0, forwardStart_)
            .withDiscountingTermStructure(discountRelinkableHandle_)
            .withEffectiveDate(startDate_)
            .withTerminationDate(endDate_)
            .withTelescopicValueDates(telescopicValueDates_)
            .withPaymentLag(paymentLag_)
            .withPaymentAdjustment(paymentConvention_)
            .withPaymentFrequency(paymentFrequency_)
            .withPaymentCalendar(paymentCalendar_)
            .withAveragingMethod(averagingMethod_)
            .withLookbackDays(lookbackDays_)
            .withLockoutDays(lockoutDays_)
            .withRule(rule_)
            .withConvention(convention_)
            .withTerminationDateConvention(convention_)
            .withObservationShift(applyObservationShift_);
        if (endOfMonth_) {
            tmp.withEndOfMonth(*endOfMonth_);
        }
        if (fixedPaymentFrequency_) {
            tmp.withFixedLegPaymentFrequency(*fixedPaymentFrequency_);
        }
        if (!fixedCalendar_.empty()) {
            tmp.withFixedLegCalendar(fixedCalendar_);
        }
        if (!overnightCalendar_.empty()) {
            tmp.withOvernightLegCalendar(overnightCalendar_);
        }
        // only set settlementDays when no explicit start date, to avoid conflict
        if (startDate_ == Date() && settlementDays_ != Null<Natural>())
            tmp.withSettlementDays(settlementDays_);
        swap_ = tmp;

        if (pricer_)
            setCouponPricer(swap_->overnightLeg(), pricer_);

        simplifyNotificationGraph(*swap_, true);

        earliestDate_ = swap_->startDate();
        maturityDate_ = swap_->maturityDate();

        Date lastPaymentDate = std::max(swap_->overnightLeg().back()->date(),
                                        swap_->fixedLeg().back()->date());
        Date lastFixingDate =
            ext::dynamic_pointer_cast<OvernightIndexedCoupon>(swap_->overnightLeg().back())->fixingDate();
        Date fixingEndDate =
            overnightIndex_->maturityDate(overnightIndex_->valueDate(lastFixingDate));
        latestRelevantDate_ = latestDate_ = std::max({maturityDate_, lastPaymentDate, fixingEndDate});

        switch (pillarChoice_) {
          case Pillar::MaturityDate:
            pillarDate_ = maturityDate_;
            break;
          case Pillar::LastRelevantDate:
            pillarDate_ = latestRelevantDate_;
            break;
          case Pillar::CustomDate:
            // pillarDate_ already assigned at construction time
            QL_REQUIRE(pillarDate_ >= earliestDate_,
                       "pillar date (" << pillarDate_ << ") must be later "
                       "than or equal to the instrument's earliest date (" <<
                       earliestDate_ << ")");
            QL_REQUIRE(pillarDate_ <= latestRelevantDate_,
                       "pillar date (" << pillarDate_ << ") must be before "
                       "or equal to the instrument's latest relevant date (" <<
                       latestRelevantDate_ << ")");
            break;
          default:
            QL_FAIL("unknown Pillar::Choice(" << Integer(pillarChoice_) << ")");
        }
    }

    void OISRateHelper::setTermStructure(YieldTermStructure* t) {
        // do not set the relinkable handle as an observer -
        // force recalculation when needed
        bool observer = false;

        ext::shared_ptr<YieldTermStructure> temp(t, null_deleter());
        termStructureHandle_.linkTo(temp, observer);

        if (discountHandle_.empty())
            discountRelinkableHandle_.linkTo(temp, observer);
        else
            discountRelinkableHandle_.linkTo(*discountHandle_, observer);

        RelativeDateRateHelper::setTermStructure(t);
    }

    Real OISRateHelper::impliedQuote() const {
        QL_REQUIRE(termStructure_ != nullptr, "term structure not set");
        // we didn't register as observers - force calculation
        swap_->deepUpdate();
        // weak implementation... to be improved
        static const Spread basisPoint = 1.0e-4;
        Real floatingLegNPV = swap_->overnightLegNPV();
        Spread spread = overnightSpread_.empty() ? 0.0 : overnightSpread_->value();
        Real spreadNPV = swap_->overnightLegBPS()/basisPoint*spread;
        Real totNPV = - (floatingLegNPV+spreadNPV);
        Real result = totNPV/(swap_->fixedLegBPS()/basisPoint);
        return result;
    }

    std::vector<std::pair<Time, Real>> OISRateHelper::impliedQuoteSensitivities() const {
        if (termStructure_ == nullptr || discountRelinkableHandle_.empty())
            return {};
        // fall back to numerical differentiation for the features that
        // the analytical formulas below don't cover
        if (averagingMethod_ != RateAveraging::Compound || lockoutDays_ != 0 ||
            applyObservationShift_ || lookbackDays_ != Null<Natural>() ||
            pricer_ != nullptr)
            return {};

        // Same structure as SwapRateHelper::impliedQuoteSensitivities();
        // the difference is in the sensitivity of the overnight coupon
        // rates, whose compound factor is (mirroring the corresponding
        // computation in CompoundingOvernightIndexedCouponPricer)
        //
        //   C = C_past * [boundary growth factors] * P(v_s)/P(v_e)
        //
        // with at most one non-telescoping growth factor at each end of
        // the value-date schedule.
        const bool sameCurve = discountHandle_.empty();
        const YieldTermStructure& discountCurve = **discountRelinkableHandle_;
        Date settlement = discountCurve.referenceDate();
        Date today = Settings::instance().evaluationDate();
        Spread s = overnightSpread_.empty() ? 0.0 : overnightSpread_->value();

        Real annuity = 0.0;
        std::vector<std::pair<Date, Real>> fixedData; // (payment date, N*tau)
        for (const auto& cf : swap_->fixedLeg()) {
            if (cf->hasOccurred(settlement))
                continue;
            auto cpn = ext::dynamic_pointer_cast<FixedRateCoupon>(cf);
            if (cpn == nullptr)
                return {};
            Real ntau = cpn->nominal()*cpn->accrualPeriod();
            annuity += ntau*discountCurve.discount(cpn->date());
            fixedData.emplace_back(cpn->date(), ntau);
        }
        if (annuity == 0.0)
            return {};

        struct OvernightData {
            Date payDate;
            Real ntau, rate, amountScale;
            DiscountFactor discount;
            // log-derivatives of the compound factor w.r.t. the
            // forecast discount factors at the given dates
            std::vector<std::pair<Date, Real>> dlogCompound;
        };
        std::vector<OvernightData> overnightData;
        Real floatingNPV = 0.0;
        const DayCounter& dc = overnightIndex_->dayCounter();
        for (const auto& cf : swap_->overnightLeg()) {
            if (cf->hasOccurred(settlement))
                continue;
            auto cpn = ext::dynamic_pointer_cast<OvernightIndexedCoupon>(cf);
            if (cpn == nullptr || cpn->compoundSpreadDaily())
                return {};

            const auto& fixingDates = cpn->fixingDates();
            const auto& valueDates = cpn->valueDates();
            const auto& interestDates = cpn->interestDates();
            const auto& dt = cpn->dt();
            Size n = fixingDates.size();

            // first fixing to be forecast
            Size i0 = 0;
            while (i0 < n && fixingDates[i0] < today)
                ++i0;
            if (i0 < n && fixingDates[i0] == today &&
                overnightIndex_->hasHistoricalFixing(fixingDates[i0]))
                ++i0;

            OvernightData d;
            if (i0 < n) {
                auto addGrowthFactor = [&](Size i) {
                    // g = 1 + f*dt[i], f = (P(v1)/P(v2)-1)/tauIndex
                    Real f = overnightIndex_->fixing(fixingDates[i]);
                    Real g = 1.0 + f*dt[i];
                    const Date& v1 = valueDates[i];
                    const Date& v2 = valueDates[i+1];
                    Time tauIndex = dc.yearFraction(v1, v2);
                    DiscountFactor P1 = termStructure_->discount(v1);
                    DiscountFactor P2 = termStructure_->discount(v2);
                    d.dlogCompound.emplace_back(v1, dt[i]/(tauIndex*P2*g));
                    d.dlogCompound.emplace_back(v2, -dt[i]*P1/(tauIndex*P2*P2*g));
                };
                // possible front and back non-telescoping factors, when
                // the accrual start or end falls on a fixing holiday
                Size start = (i0 == 0 && valueDates.front() < interestDates.front())
                             ? Size(1) : i0;
                Size end = n - (valueDates[n] <= interestDates[n] ? 0 : 1);
                if (start < end) {
                    for (Size i = i0; i < start; ++i)
                        addGrowthFactor(i);
                    DiscountFactor Ps = termStructure_->discount(valueDates[start]);
                    DiscountFactor Pe = termStructure_->discount(valueDates[end]);
                    d.dlogCompound.emplace_back(valueDates[start], 1.0/Ps);
                    d.dlogCompound.emplace_back(valueDates[end], -1.0/Pe);
                    for (Size i = end; i < n; ++i)
                        addGrowthFactor(i);
                } else {
                    for (Size i = i0; i < n; ++i)
                        addGrowthFactor(i);
                }
            }

            d.payDate = cpn->date();
            d.ntau = cpn->nominal()*cpn->accrualPeriod();
            d.rate = cpn->rate();
            d.discount = discountCurve.discount(d.payDate);
            // amount = N * accrualPeriod * gearing * (compound-1)/tau
            Real gearing = cpn->gearing();
            Time tau = dc.yearFraction(interestDates.front(), interestDates.back());
            Real compound = (d.rate - cpn->spread())*tau/gearing + 1.0;
            d.amountScale = d.ntau*gearing*compound/tau;

            floatingNPV += d.ntau*(d.rate + s)*d.discount;
            overnightData.push_back(std::move(d));
        }
        Real fairRate = floatingNPV/annuity;

        std::vector<std::pair<Time, Real>> result;
        if (sameCurve) {
            // sensitivities to the discount factors at the payment dates
            for (const auto& [payDate, ntau] : fixedData)
                result.emplace_back(termStructure_->timeFromReference(payDate),
                                    -fairRate*ntau/annuity);
            for (const auto& d : overnightData)
                result.emplace_back(termStructure_->timeFromReference(d.payDate),
                                    d.ntau*(d.rate + s)/annuity);
        }
        // sensitivities to the forecast discount factors
        for (const auto& d : overnightData)
            for (const auto& [date, w] : d.dlogCompound)
                result.emplace_back(termStructure_->timeFromReference(date),
                                    d.discount*d.amountScale*w/annuity);
        return result;
    }

    void OISRateHelper::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<OISRateHelper>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            RateHelper::accept(v);
    }

}
