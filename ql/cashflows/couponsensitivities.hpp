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
#include <ql/time/date.hpp>
#include <ql/types.hpp>
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

        /*! Sensitivities of the fair fixed rate of a fixed-vs-floating
            swap (plus a spread quoted on the floating leg) to the
            discount factors of the curve being bootstrapped; the
            floating-leg forwards are assumed to be forecast on that
            curve, while the discount curve may or may not be the same
            curve.  An empty vector means that some coupon is not
            supported by the analytical formulas.
        */
        std::vector<std::pair<Time, Real>>
        fairRateSensitivities(const Leg& fixedLeg,
                              const Leg& floatingLeg,
                              Spread helperSpread,
                              const YieldTermStructure* bootstrappedCurve,
                              const YieldTermStructure& discountCurve,
                              bool discountOnBootstrappedCurve);

        /*! Sensitivities of the fair basis of a floating-vs-floating
            swap paying baseLeg + basis and receiving otherLeg, i.e., of
            \f$ b = (O - B)/A \f$ with O and B the leg NPVs and A the
            base-leg annuity, to the discount factors of the curve being
            bootstrapped.  The flags tell which parts of the pricing use
            that curve.  An empty vector means that some coupon is not
            supported by the analytical formulas.
        */
        std::vector<std::pair<Time, Real>>
        fairBasisSensitivities(const Leg& baseLeg,
                               const Leg& otherLeg,
                               const YieldTermStructure* bootstrappedCurve,
                               const YieldTermStructure& discountCurve,
                               bool discountOnBootstrappedCurve,
                               bool baseLegForecastsOnBootstrappedCurve,
                               bool otherLegForecastsOnBootstrappedCurve);

    }

}

#endif
