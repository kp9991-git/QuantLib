/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Quantlib contributors

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
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/settings.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantLib {

    namespace detail {

        CouponSensitivityAnalysis analyzeCoupon(const ext::shared_ptr<CashFlow>& cf,
                                                bool withSensitivities) {
            CouponSensitivityAnalysis a;

            if (!withSensitivities) {
                // when only the value is needed, any coupon will do
                if (auto cpn = ext::dynamic_pointer_cast<Coupon>(cf)) {
                    a.supported = true;
                    a.ntau = cpn->nominal()*cpn->accrualPeriod();
                    a.amount = cpn->amount();
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

                // The compound factor is (mirroring the corresponding
                // computation in CompoundingOvernightIndexedCouponPricer)
                //
                //   C = C_past * [boundary growth factors] * P(v_s)/P(v_e)
                //
                // with at most one non-telescoping growth factor at each
                // end of the value-date schedule.
                auto index = ext::dynamic_pointer_cast<OvernightIndex>(overnight->index());
                const Handle<YieldTermStructure>& curve = index->forwardingTermStructure();
                const DayCounter& dc = index->dayCounter();
                Date today = Settings::instance().evaluationDate();

                const auto& fixingDates = overnight->fixingDates();
                const auto& valueDates = overnight->valueDates();
                const auto& interestDates = overnight->interestDates();
                const auto& dt = overnight->dt();
                Size n = fixingDates.size();

                // first fixing to be forecast
                Size i0 = 0;
                while (i0 < n && fixingDates[i0] < today)
                    ++i0;
                if (i0 < n && fixingDates[i0] == today &&
                    index->hasHistoricalFixing(fixingDates[i0]))
                    ++i0;
                if (i0 == n)  // fully fixed
                    return a;

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
                // possible front and back non-telescoping factors, when
                // the accrual start or end falls on a fixing holiday
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

        std::vector<std::pair<Time, Real>>
        fairRateSensitivities(const Leg& fixedLeg,
                              const Leg& floatingLeg,
                              Spread helperSpread,
                              const YieldTermStructure* bootstrappedCurve,
                              const YieldTermStructure& discountCurve,
                              bool discountOnBootstrappedCurve) {
            Date settlement = discountCurve.referenceDate();

            // fair rate = F/A, with F the floating-leg NPV including the
            // helper spread and A the fixed-leg annuity
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

            struct FloatingData {
                Date payDate;
                Real ntau, amount;
                DiscountFactor discount;
                std::vector<std::pair<Date, Real>> sensitivities;
            };
            std::vector<FloatingData> floatingData;
            Real floatingNPV = 0.0;
            for (const auto& cf : floatingLeg) {
                if (cf->hasOccurred(settlement))
                    continue;
                auto a = analyzeCoupon(cf, true);
                if (!a.supported)
                    return {};
                FloatingData d = {cf->date(), a.ntau, a.amount,
                                  discountCurve.discount(cf->date()),
                                  std::move(a.amountSensitivities)};
                floatingNPV += (d.amount + helperSpread*d.ntau)*d.discount;
                floatingData.push_back(std::move(d));
            }
            Real fairRate = floatingNPV/annuity;

            std::vector<std::pair<Time, Real>> result;
            if (discountOnBootstrappedCurve) {
                // sensitivities to the discount factors at the payment dates
                for (const auto& [payDate, ntau] : fixedData)
                    result.emplace_back(bootstrappedCurve->timeFromReference(payDate),
                                        -fairRate*ntau/annuity);
                for (const auto& d : floatingData)
                    result.emplace_back(bootstrappedCurve->timeFromReference(d.payDate),
                                        (d.amount + helperSpread*d.ntau)/annuity);
            }
            // sensitivities to the forecast discount factors
            for (const auto& d : floatingData)
                for (const auto& [date, w] : d.sensitivities)
                    result.emplace_back(bootstrappedCurve->timeFromReference(date),
                                        d.discount*w/annuity);
            return result;
        }

        std::vector<std::pair<Time, Real>>
        fairBasisSensitivities(const Leg& baseLeg,
                               const Leg& otherLeg,
                               const YieldTermStructure* bootstrappedCurve,
                               const YieldTermStructure& discountCurve,
                               bool discountOnBootstrappedCurve,
                               bool baseLegForecastsOnBootstrappedCurve,
                               bool otherLegForecastsOnBootstrappedCurve) {
            Date settlement = discountCurve.referenceDate();

            struct LegData {
                Date payDate;
                Real ntau, amount;
                DiscountFactor discount;
                std::vector<std::pair<Date, Real>> sensitivities;
            };
            auto analyzeLeg = [&](const Leg& leg, bool withSensitivities,
                                  std::vector<LegData>& data) {
                for (const auto& cf : leg) {
                    if (cf->hasOccurred(settlement))
                        continue;
                    auto a = analyzeCoupon(cf, withSensitivities);
                    if (!a.supported)
                        return false;
                    data.push_back({cf->date(), a.ntau, a.amount,
                                    discountCurve.discount(cf->date()),
                                    std::move(a.amountSensitivities)});
                }
                return true;
            };

            std::vector<LegData> baseData, otherData;
            if (!analyzeLeg(baseLeg, baseLegForecastsOnBootstrappedCurve, baseData) ||
                !analyzeLeg(otherLeg, otherLegForecastsOnBootstrappedCurve, otherData))
                return {};

            // fair basis = (O - B)/A, with B and O the leg NPVs and A the
            // base-leg annuity
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

            std::vector<std::pair<Time, Real>> result;
            if (discountOnBootstrappedCurve) {
                // sensitivities to the discount factors at the payment dates
                for (const auto& d : baseData)
                    result.emplace_back(bootstrappedCurve->timeFromReference(d.payDate),
                                        (-d.amount - fairBasis*d.ntau)/annuity);
                for (const auto& d : otherData)
                    result.emplace_back(bootstrappedCurve->timeFromReference(d.payDate),
                                        d.amount/annuity);
            }
            // sensitivities to the forecast discount factors
            if (baseLegForecastsOnBootstrappedCurve)
                for (const auto& d : baseData)
                    for (const auto& [date, w] : d.sensitivities)
                        result.emplace_back(bootstrappedCurve->timeFromReference(date),
                                            -d.discount*w/annuity);
            if (otherLegForecastsOnBootstrappedCurve)
                for (const auto& d : otherData)
                    for (const auto& [date, w] : d.sensitivities)
                        result.emplace_back(bootstrappedCurve->timeFromReference(date),
                                            d.discount*w/annuity);
            return result;
        }

    }

}
