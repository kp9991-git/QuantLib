/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

#include <ql/cashflows/couponsensitivities.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/floatingratecoupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/settings.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantLib {

    namespace detail {

        namespace {

            const YieldTermStructure* couponForecastCurve(const ext::shared_ptr<CashFlow>& cf) {
                auto floating = ext::dynamic_pointer_cast<FloatingRateCoupon>(cf);
                if (floating == nullptr)
                    return nullptr;
                auto index = ext::dynamic_pointer_cast<IborIndex>(floating->index());
                if (index == nullptr)
                    return nullptr;
                const Handle<YieldTermStructure>& h = index->forwardingTermStructure();
                return h.empty() ? nullptr : &**h;
            }

        }

        CouponSensitivityAnalysis analyzeCoupon(const ext::shared_ptr<CashFlow>& cf,
                                                bool withSensitivities) {
            CouponSensitivityAnalysis a;

            if (!withSensitivities) {
                // amount-only mode supports any coupon
                if (auto cpn = ext::dynamic_pointer_cast<Coupon>(cf)) {
                    a.supported = true;
                    a.ntau = cpn->nominal()*cpn->accrualPeriod();
                    a.amount = cpn->amount();
                    a.forecastCurve = couponForecastCurve(cf);
                }
                return a;
            }

            if (auto fixed = ext::dynamic_pointer_cast<FixedRateCoupon>(cf)) {
                a.supported = true;
                a.ntau = fixed->nominal()*fixed->accrualPeriod();
                a.amount = fixed->amount();
                return a;
            }

            if (auto ibor = ext::dynamic_pointer_cast<IborCoupon>(cf)) {
                if (ibor->isInArrears())
                    return a;
                a.supported = true;
                a.ntau = ibor->nominal()*ibor->accrualPeriod();
                Real gearing = ibor->gearing();
                a.amount = a.ntau*(gearing*ibor->indexFixing() + ibor->spread());
                if (withSensitivities && !ibor->hasFixed()) {
                    // f = (P(v)/P(e)-1)/tau on the index's forecast curve
                    const Handle<YieldTermStructure>& curve =
                        ibor->iborIndex()->forwardingTermStructure();
                    a.forecastCurve = &**curve;
                    const Date& v = ibor->fixingValueDate();
                    const Date& e = ibor->fixingEndDate();
                    Time tau = ibor->spanningTime();
                    DiscountFactor Pv = curve->discount(v);
                    DiscountFactor Pe = curve->discount(e);
                    Real k = gearing*a.ntau/tau;
                    a.amountSensitivities.emplace_back(v, k/Pe);
                    a.amountSensitivities.emplace_back(e, -k*Pv/(Pe*Pe));
                }
                return a;
            }

            if (auto overnight = ext::dynamic_pointer_cast<OvernightIndexedCoupon>(cf)) {
                if (overnight->averagingMethod() != RateAveraging::Compound ||
                    overnight->compoundSpreadDaily() ||
                    overnight->lockoutDays() != 0 ||
                    overnight->applyObservationShift() ||
                    !overnight->canApplyTelescopicFormula())
                    return a;
                a.supported = true;
                a.ntau = overnight->nominal()*overnight->accrualPeriod();
                Real rate = overnight->rate();
                a.amount = a.ntau*rate;
                if (!withSensitivities)
                    return a;

                // Same factor as CompoundingOvernightIndexedCouponPricer
                //   C = C_past * [boundary growth factors] * P(v_s)/P(v_e)
                // with at most one boundary factor at each end
                auto index = ext::dynamic_pointer_cast<OvernightIndex>(overnight->index());
                const Handle<YieldTermStructure>& curve = index->forwardingTermStructure();
                const DayCounter& dc = index->dayCounter();
                Date today = Settings::instance().evaluationDate();

                const auto& fixingDates = overnight->fixingDates();
                const auto& valueDates = overnight->valueDates();
                const auto& interestDates = overnight->interestDates();
                const auto& dt = overnight->dt();
                Size n = fixingDates.size();

                // first forecast fixing
                Size i0 = 0;
                while (i0 < n && fixingDates[i0] < today)
                    ++i0;
                if (i0 < n && fixingDates[i0] == today &&
                    index->hasHistoricalFixing(fixingDates[i0]))
                    ++i0;
                if (i0 == n)  // fully fixed
                    return a;
                a.forecastCurve = &**curve;

                Real gearing = overnight->gearing();
                Time tau = dc.yearFraction(interestDates.front(), interestDates.back());
                Real compound = (rate - overnight->spread())*tau/gearing + 1.0;
                // amount = ntau * (gearing*(C-1)/tau + spread)
                Real amountScale = a.ntau*gearing*compound/tau;

                auto addGrowthFactor = [&](Size i) {
                    // g = 1 + f*dt[i], f = (P(v1)/P(v2)-1)/tauIndex
                    Real f = index->fixing(fixingDates[i]);
                    Real g = 1.0 + f*dt[i];
                    const Date& v1 = valueDates[i];
                    const Date& v2 = valueDates[i+1];
                    Time tauIndex = dc.yearFraction(v1, v2);
                    DiscountFactor P1 = curve->discount(v1);
                    DiscountFactor P2 = curve->discount(v2);
                    a.amountSensitivities.emplace_back(
                        v1, amountScale*dt[i]/(tauIndex*P2*g));
                    a.amountSensitivities.emplace_back(
                        v2, -amountScale*dt[i]*P1/(tauIndex*P2*P2*g));
                };
                // boundary factors caused by fixing holidays
                Size start = (i0 == 0 && valueDates.front() < interestDates.front())
                             ? Size(1) : i0;
                Size end = n - (valueDates[n] <= interestDates[n] ? 0 : 1);
                if (start < end) {
                    for (Size i = i0; i < start; ++i)
                        addGrowthFactor(i);
                    DiscountFactor Ps = curve->discount(valueDates[start]);
                    DiscountFactor Pe = curve->discount(valueDates[end]);
                    a.amountSensitivities.emplace_back(valueDates[start],
                                                       amountScale/Ps);
                    a.amountSensitivities.emplace_back(valueDates[end],
                                                       -amountScale/Pe);
                    for (Size i = end; i < n; ++i)
                        addGrowthFactor(i);
                } else {
                    for (Size i = i0; i < n; ++i)
                        addGrowthFactor(i);
                }
                return a;
            }

            return a;
        }

        bool analyzeFloatingLeg(const Leg& leg,
                                const Date& settlement,
                                const YieldTermStructure& discountCurve,
                                std::vector<FloatingFlowData>& data,
                                QuoteSensitivities& result,
                                std::optional<bool> includeSettlementDateFlows) {
            for (const auto& cf : leg) {
                if (cf->hasOccurred(settlement, includeSettlementDateFlows))
                    continue;
                auto a = analyzeCoupon(cf, true);
                if (!a.supported) {
                    a = analyzeCoupon(cf, false);
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
                data.push_back({cf->date(), a.ntau, a.amount,
                                discountCurve.discount(cf->date()),
                                a.forecastCurve,
                                std::move(a.amountSensitivities)});
            }
            return true;
        }

        QuoteSensitivities
        fairRateSensitivities(const Leg& fixedLeg,
                              const Leg& floatingLeg,
                              Spread helperSpread,
                              const YieldTermStructure& discountCurve) {
            QuoteSensitivities result;
            Date settlement = discountCurve.referenceDate();

            // fair rate = floating-leg NPV / fixed-leg annuity
            Real annuity = 0.0;
            std::vector<std::pair<Date, Real>> fixedData; // (payment date, N*tau)
            for (const auto& cf : fixedLeg) {
                if (cf->hasOccurred(settlement))
                    continue;
                auto a = analyzeCoupon(cf, false);
                if (!a.supported)
                    return {};
                annuity += a.ntau*discountCurve.discount(cf->date());
                fixedData.emplace_back(cf->date(), a.ntau);
            }
            if (annuity == 0.0)
                return {};

            std::vector<FloatingFlowData> floatingData;
            if (!analyzeFloatingLeg(floatingLeg, settlement, discountCurve,
                                    floatingData, result))
                return {};
            Real floatingNPV = 0.0;
            for (const auto& d : floatingData)
                floatingNPV += (d.amount + helperSpread*d.ntau)*d.discount;
            Real fairRate = floatingNPV/annuity;

            // discount sensitivities
            auto& discountBucket =
                result.sensitivities[static_cast<const TermStructure*>(&discountCurve)];
            for (const auto& [payDate, ntau] : fixedData)
                discountBucket.emplace_back(payDate, -fairRate*ntau/annuity);
            for (const auto& d : floatingData)
                discountBucket.emplace_back(d.payDate,
                                            (d.amount + helperSpread*d.ntau)/annuity);
            // forecast sensitivities
            for (const auto& d : floatingData) {
                if (d.sensitivities.empty())
                    continue;
                auto& bucket = result.sensitivities[
                    static_cast<const TermStructure*>(d.forecastCurve)];
                for (const auto& [date, w] : d.sensitivities)
                    bucket.emplace_back(date, d.discount*w/annuity);
            }
            result.available = true;
            return result;
        }

        QuoteSensitivities
        fairBasisSensitivities(const Leg& baseLeg,
                               const Leg& otherLeg,
                               const YieldTermStructure& discountCurve) {
            QuoteSensitivities result;
            Date settlement = discountCurve.referenceDate();

            std::vector<FloatingFlowData> baseData, otherData;
            if (!analyzeFloatingLeg(baseLeg, settlement, discountCurve,
                                    baseData, result) ||
                !analyzeFloatingLeg(otherLeg, settlement, discountCurve,
                                    otherData, result))
                return {};

            // fair basis = (other NPV - base NPV) / base annuity
            Real annuity = 0.0, baseNPV = 0.0, otherNPV = 0.0;
            for (const auto& d : baseData) {
                annuity += d.ntau*d.discount;
                baseNPV += d.amount*d.discount;
            }
            for (const auto& d : otherData)
                otherNPV += d.amount*d.discount;
            if (annuity == 0.0)
                return {};
            Real fairBasis = (otherNPV - baseNPV)/annuity;

            // discount sensitivities
            auto& discountBucket =
                result.sensitivities[static_cast<const TermStructure*>(&discountCurve)];
            for (const auto& d : baseData)
                discountBucket.emplace_back(d.payDate,
                                            (-d.amount - fairBasis*d.ntau)/annuity);
            for (const auto& d : otherData)
                discountBucket.emplace_back(d.payDate, d.amount/annuity);
            // forecast sensitivities
            for (const auto& d : baseData) {
                if (d.sensitivities.empty())
                    continue;
                auto& bucket = result.sensitivities[
                    static_cast<const TermStructure*>(d.forecastCurve)];
                for (const auto& [date, w] : d.sensitivities)
                    bucket.emplace_back(date, -d.discount*w/annuity);
            }
            for (const auto& d : otherData) {
                if (d.sensitivities.empty())
                    continue;
                auto& bucket = result.sensitivities[
                    static_cast<const TermStructure*>(d.forecastCurve)];
                for (const auto& [date, w] : d.sensitivities)
                    bucket.emplace_back(date, d.discount*w/annuity);
            }
            result.available = true;
            return result;
        }

    }

}
