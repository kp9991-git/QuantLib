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

/*! \file couponsensitivities.hpp
    \brief analytical sensitivities of coupon amounts and fair rates,
           used by the rate helpers to provide analytical Jacobians
*/

#ifndef quantlib_coupon_sensitivities_hpp
#define quantlib_coupon_sensitivities_hpp

#include <ql/cashflow.hpp>
#include <ql/termstructures/bootstraphelper.hpp>
#include <ql/time/date.hpp>
#include <ql/types.hpp>
#include <optional>
#include <utility>
#include <vector>

namespace QuantLib {

    class YieldTermStructure;

    namespace detail {

        //! analytical sensitivity analysis of a single coupon
        struct CouponSensitivityAnalysis {
            bool supported = false;
            //! nominal times accrual period
            Real ntau = 0.0;
            //! coupon amount, consistent with the pricing engines
            Real amount = 0.0;
            //! the coupon's forecast curve, when the coupon has one
            //! that could be identified
            const YieldTermStructure* forecastCurve = nullptr;
            //! derivative of the amount with respect to the discount
            //! factors of the coupon's forecast curve at the given dates
            std::vector<std::pair<Date, Real>> amountSensitivities;
        };

        /*! Analyzes fixed-rate coupons, vanilla ibor coupons, and
            compounded overnight coupons.  Other cash-flow types, or
            coupons with features not covered by the analytical formulas
            (in-arrears fixings, daily compounded spreads, lookback,
            lockout, observation shift, arithmetic averaging) are
            flagged as not supported.
        */
        CouponSensitivityAnalysis analyzeCoupon(const ext::shared_ptr<CashFlow>& cf,
                                                bool withSensitivities);

        //! analysis of a single floating coupon, discounted on a given curve
        struct FloatingFlowData {
            Date payDate;
            Real ntau, amount;
            DiscountFactor discount;
            const YieldTermStructure* forecastCurve;
            std::vector<std::pair<Date, Real>> sensitivities;
        };

        /*! Analyzes the coupons of a floating leg, trying the full
            analytical formulas first and falling back to pricing by
            amount() when they don't apply; in the latter case the
            coupon's forecast curve is added to the incomplete set of
            the given result.  Returns false when the leg can't be
            analyzed at all.
        */
        bool analyzeFloatingLeg(const Leg& leg,
                                const Date& settlement,
                                const YieldTermStructure& discountCurve,
                                std::vector<FloatingFlowData>& data,
                                QuoteSensitivities& result,
                                std::optional<bool> includeSettlementDateFlows = std::nullopt);

        /*! Sensitivities of the fair fixed rate of a fixed-vs-floating
            swap (plus a spread quoted on the floating leg) to the
            discount factors of every curve entering the pricing: the
            given discount curve and the forecast curves of the
            floating-leg coupons.  Floating coupons not supported by the
            analytical formulas are priced by amount() and their
            forecast curves are reported as incomplete; if that fails
            too, the whole result is flagged as unavailable.
        */
        QuoteSensitivities
        fairRateSensitivities(const Leg& fixedLeg,
                              const Leg& floatingLeg,
                              Spread helperSpread,
                              const YieldTermStructure& discountCurve);

        /*! Sensitivities of the fair basis of a floating-vs-floating
            swap paying baseLeg + basis and receiving otherLeg, i.e., of
            \f$ b = (O - B)/A \f$ with O and B the leg NPVs and A the
            base-leg annuity, to the discount factors of every curve
            entering the pricing: the given discount curve and the
            forecast curves of the coupons on both legs.  Coupons not
            supported by the analytical formulas are priced by amount()
            and their forecast curves are reported as incomplete; if
            that fails too, the whole result is flagged as unavailable.
        */
        QuoteSensitivities
        fairBasisSensitivities(const Leg& baseLeg,
                               const Leg& otherLeg,
                               const YieldTermStructure& discountCurve);

    }

}

#endif
