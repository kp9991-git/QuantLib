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

/*! \file quotesensitivitycalculator.hpp
    \brief building blocks for analytical implied-quote sensitivities
*/

#ifndef quantlib_quote_sensitivity_calculator_hpp
#define quantlib_quote_sensitivity_calculator_hpp

#include <ql/cashflow.hpp>
#include <ql/termstructures/bootstraphelper.hpp>
#include <ql/time/date.hpp>
#include <ql/types.hpp>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace QuantLib {

    class YieldTermStructure;

    namespace detail {

        //! sensitivity to one curve value
        struct CurvePointSensitivity {
            const TermStructure* curve;
            Date date;
            Real derivative;
        };
        using CurvePointSensitivities = std::vector<CurvePointSensitivity>;

        //! Adds sensitivities of \f$ s(P(d_1)/P(d_2)-1)/\tau \f$
        void addSimpleForwardSensitivities(QuoteSensitivities& result,
                                           const YieldTermStructure& curve,
                                           const Date& d1,
                                           const Date& d2,
                                           Time tau,
                                           Real scale = 1.0);

        //! coupon value and forecast sensitivities
        struct CouponSensitivityAnalysis {
            bool supported = false;
            //! nominal times accrual
            Real ntau = 0.0;
            //! coupon amount
            Real amount = 0.0;
            //! single forecast curve used by unsupported-coupon fallback
            const YieldTermStructure* forecastCurve = nullptr;
            //! amount derivatives with respect to forecast discount factors
            CurvePointSensitivities amountSensitivities;
        };

        /*! Analyzes fixed, Ibor (including weighted-index stubs), and
            compounded overnight coupons. Unsupported cases return
            supported=false. */
        CouponSensitivityAnalysis analyzeCoupon(const ext::shared_ptr<CashFlow>& cf,
                                                bool withSensitivities);

        /*! Retries value-only when analytic sensitivities are unavailable
            and marks the forecast curve incomplete. */
        bool analyzeCouponWithFallback(const ext::shared_ptr<CashFlow>& cf,
                                       CouponSensitivityAnalysis& analysis,
                                       QuoteSensitivities& result);

        //! undiscounted cash-flow data
        struct FlowSensitivityData {
            Date payDate;
            Real amount = 0.0;
            //! nominal times accrual period (zero for non-coupon flows)
            Real ntau = 0.0;
            //! amount and ntau derivatives
            CurvePointSensitivities amountSensitivities, ntauSensitivities;
        };

        //! discounted leg values and their curve sensitivities
        struct LegSensitivityAnalysis {
            const YieldTermStructure* discountCurve = nullptr;
            Real npv = 0.0;
            Real annuity = 0.0;
            CurvePointSensitivities npvSensitivities, annuitySensitivities;

            //! discounts and adds the flow
            void addFlow(const FlowSensitivityData& flow);
            //! adds a flow with curve-independent values
            void addFlow(const Date& payDate, Real amount, Real ntau = 0.0);
        };

        //! result from a custom cash-flow analyzer
        enum class FlowHandling {
            NotApplicable, //!< use the default analysis
            Analyzed,      //!< data supplied
            Unsupported    //!< cannot analyze the leg
        };
        using FlowAnalyzer =
            std::function<FlowHandling(const ext::shared_ptr<CashFlow>&,
                                       FlowSensitivityData&)>;

        /*! Adds future cash flows to the analysis, trying the custom analyzer
            first. Returns false if a flow cannot be analyzed. */
        bool analyzeLeg(const Leg& leg,
                        const YieldTermStructure& discountCurve,
                        QuoteSensitivities& result,
                        LegSensitivityAnalysis& analysis,
                        std::optional<bool> includeSettlementDateFlows = std::nullopt,
                        const FlowAnalyzer& customFlows = {});

        //! weighted NPV and annuity contribution
        struct LegTerm {
            const LegSensitivityAnalysis* leg;
            Real amountWeight = 0.0;
            Real ntauWeight = 0.0;
        };

        //! \f$ Q=N/D \f$ from weighted legs and optional extra derivatives
        struct QuotientSensitivitySpec {
            std::vector<LegTerm> numerator, denominator;
            CurvePointSensitivities numeratorExtra, denominatorExtra;
        };

        //! Adds \f$ dQ=dN/D-Q\,dD/D \f$ or returns false when \f$ D=0 \f$
        bool addQuotientSensitivities(QuoteSensitivities& result,
                                      const QuotientSensitivitySpec& spec);

        //! sensitivities of a fixed-vs-floating fair rate
        QuoteSensitivities fairRateSensitivities(
            const Leg& fixedLeg,
                              const Leg& floatingLeg,
                              Spread helperSpread,
                              const YieldTermStructure& discountCurve);

        //! sensitivities of floating-vs-floating fair basis \f$ b=(O-B)/A \f$
        QuoteSensitivities fairBasisSensitivities(
            const Leg& baseLeg,
                               const Leg& otherLeg,
                               const YieldTermStructure& discountCurve,
                               std::optional<bool> includeSettlementDateFlows = std::nullopt);

    }

}

#endif
