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

/*! \file couponsensitivities.hpp
    \brief coupon and fair-rate sensitivities for analytical Jacobians
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

        //! analytical data for one coupon
        struct CouponSensitivityAnalysis {
            bool supported = false;
            //! nominal times accrual period
            Real ntau = 0.0;
            //! coupon amount
            Real amount = 0.0;
            //! identified forecast curve, if any
            const YieldTermStructure* forecastCurve = nullptr;
            //! derivative of the amount with respect to the discount
            //! factors of the coupon's forecast curve at the given dates
            std::vector<std::pair<Date, Real>> amountSensitivities;
        };

        /*! Supports fixed, vanilla ibor, and compounded overnight coupons.
            Other types and unsupported features return supported=false.
        */
        CouponSensitivityAnalysis analyzeCoupon(const ext::shared_ptr<CashFlow>& cf,
                                                bool withSensitivities);

        //! discounted data for one floating cash flow
        struct FloatingFlowData {
            Date payDate;
            Real ntau, amount;
            DiscountFactor discount;
            const YieldTermStructure* forecastCurve;
            std::vector<std::pair<Date, Real>> sensitivities;
        };

        /*! Analyzes a floating leg. Unsupported formulas fall back to
            amount() and mark the forecast curve as incomplete. Returns
            false if the cash flow cannot be priced or its curve identified.
        */
        bool analyzeFloatingLeg(const Leg& leg,
                                const Date& settlement,
                                const YieldTermStructure& discountCurve,
                                std::vector<FloatingFlowData>& data,
                                QuoteSensitivities& result,
                                std::optional<bool> includeSettlementDateFlows = std::nullopt);

        /*! Per-curve sensitivities of a fixed-vs-floating fair rate,
            including any spread on the floating leg.
        */
        QuoteSensitivities
        fairRateSensitivities(const Leg& fixedLeg,
                              const Leg& floatingLeg,
                              Spread helperSpread,
                              const YieldTermStructure& discountCurve);

        /*! Per-curve sensitivities of the fair basis
            \f$ b=(O-B)/A \f$ of a floating-vs-floating swap.
        */
        QuoteSensitivities
        fairBasisSensitivities(const Leg& baseLeg,
                               const Leg& otherLeg,
                               const YieldTermStructure& discountCurve);

    }

}

#endif
