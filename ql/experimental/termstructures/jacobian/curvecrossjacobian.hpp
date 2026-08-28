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

/*! \file curvecrossjacobian.hpp
    \brief cross-curve Jacobian construction
*/

#ifndef quantlib_experimental_curve_cross_jacobian_hpp
#define quantlib_experimental_curve_cross_jacobian_hpp

#include <ql/experimental/termstructures/jacobian/bootstrapequationjacobian.hpp>
#include <ql/experimental/termstructures/jacobian/curvechainrulecalculator.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace QuantLib {

    namespace detail {

        template <class C, class = void>
        constexpr bool hasBaseCurveHandle = false;

        template <class C>
        constexpr bool hasBaseCurveHandle<
            C, std::void_t<decltype(std::declval<const C&>().baseCurve())>> = true;

        template <class T, class Curve, class = void>
        constexpr bool hasBaseCurveSensitivityTransform = false;

        template <class T, class Curve>
        constexpr bool hasBaseCurveSensitivityTransform<
            T, Curve,
            std::void_t<decltype(T::transformBaseCurveSensitivities(
                std::declval<const Curve*>(),
                std::declval<const DatedCurveSensitivities&>(),
                std::declval<DatedCurveSensitivities&>()))>> = true;

        //! type-erased interface to a bootstrapped curve
        /*! A null curve means the interface is unavailable. */
        struct CurveJacobianNode {
            ext::shared_ptr<YieldTermStructure> curve;
            CurveId id = nullptr;
            //! term structures used to turn nodes into curve values
            CurveDependencies valueDependencies;
            std::function<void()> ensure;
            std::function<Size()> numNodes;
            std::function<std::vector<Date>()> nodeDates;
            std::function<std::vector<ext::shared_ptr<BootstrapHelper<YieldTermStructure>>>()>
                aliveHelpers;
            std::function<Matrix(std::vector<bool>*)> ownJacobian;
            std::function<bool(const std::vector<std::pair<Date, Real>>&,
                               std::vector<Real>&)> analyticRow;
            std::function<Real(Size)> nodeValue;
            std::function<void(Size, Real)> setNodeValue;
        };

        //! calibrated curves and derived wrappers in one group
        struct CurveJacobianGroup {
            std::vector<CurveJacobianNode> members;
            //! index of the curve that requested the group
            Size target = 0;
            std::set<const TermStructure*> dependents;
        };

        // bootstrap classes that expose a multi-curve group
        template <class B, class = void>
        constexpr bool hasJacobianGroup = false;

        template <class B>
        constexpr bool hasJacobianGroup<
            B, std::void_t<decltype(std::declval<const B&>().jacobianGroup())>> = true;

        // bootstrapped yield curves exposing jacobian()
        template <class Curve, class = void>
        constexpr bool supportsCurveJacobianNode = false;

        template <class Curve>
        constexpr bool supportsCurveJacobianNode<
            Curve, std::void_t<decltype(std::declval<const Curve&>().jacobian(
                std::declval<std::vector<bool>*>()))>> =
            std::is_base_of_v<YieldTermStructure, Curve>;

        template <class Curve>
        struct BootstrapJacobianAccess {
            using Traits = typename Curve::traits_type;

            static CurveJacobianNode makeNode(const ext::shared_ptr<Curve>& curve) {
                CurveJacobianNode n;
                n.curve = curve;
                n.id = static_cast<const TermStructure*>(curve.get());
                if constexpr (hasBaseCurveHandle<Curve>) {
                    if (!curve->baseCurve().empty()) {
                        const auto* baseId = static_cast<const TermStructure*>(
                            curve->baseCurve().currentLink().get());
                        if constexpr (hasBaseCurveSensitivityTransform<Traits, Curve>) {
                            n.valueDependencies.add(
                                baseId,
                                [curve](const DatedCurveSensitivities& input,
                                        DatedCurveSensitivities& output) {
                                    return Traits::transformBaseCurveSensitivities(
                                        curve.get(), input, output);
                                });
                        } else {
                            n.valueDependencies.add(baseId);
                        }
                    }
                }
                n.ensure = [curve] { curve->calculate(); };
                n.numNodes = [curve]() -> Size {
                    curve->calculate();
                    return curve->times_.size() - 1;
                };
                n.nodeDates = [curve] {
                    curve->calculate();
                    return std::vector<Date>(curve->dates_.begin() + 1,
                                             curve->dates_.end());
                };
                n.aliveHelpers = [curve] {
                    curve->calculate();
                    Date firstDate = Traits::initialDate(curve.get());
                    std::vector<ext::shared_ptr<BootstrapHelper<YieldTermStructure>>> alive;
                    for (const auto& helper : curve->instruments_)
                        if (helper->pillarDate() > firstDate)
                            alive.push_back(helper);
                    return alive;
                };
                n.ownJacobian = [curve](std::vector<bool>* analyticEquations) {
                    return curve->jacobian(analyticEquations);
                };
                n.analyticRow = [curve](const std::vector<std::pair<Date, Real>>& dateSens,
                                        std::vector<Real>& row) {
                    curve->calculate();
                    // analytical weights do not include jumps
                    if (!curve->jumpDates().empty())
                        return false;
                    std::vector<std::pair<Time, Real>> sensitivities;
                    sensitivities.reserve(dateSens.size());
                    for (const auto& [date, dQdP] : dateSens)
                        sensitivities.emplace_back(curve->timeFromReference(date), dQdP);
                    return analyticBootstrapEquationRow<Traits>(
                        curve.get(), sensitivities, curve->times_,
                        curve->interpolation_, row);
                };
                n.nodeValue = [curve](Size j) {
                    curve->calculate();
                    return curve->data_[j];
                };
                n.setNodeValue = [curve](Size j, Real v) {
                    // preserve interpolation iterators
                    curve->jacobianCache_.invalidate();
                    Traits::updateGuess(curve->data_, v, j);
                    curve->interpolation_.update();
                };
                return n;
            }
        };

        //! Fill selected rows by perturbing each free node
        template <class Value>
        void fillNumericalNodeRows(const CurveJacobianNode& n,
                                   const std::vector<Size>& rows,
                                   Matrix& J,
                                   const Value& value) {
            if (rows.empty())
                return;
            std::vector<Real> up(rows.size());
            for (Size j = 1; j <= J.columns(); ++j) {
                Real v = n.nodeValue(j);
                Real h = 1.0e-6 * std::max(std::abs(v), 0.01);

                struct RestoreNode {  // NOLINT(cppcoreguidelines-special-member-functions)
                    const CurveJacobianNode& node;
                    Size j;
                    Real value;
                    ~RestoreNode() { node.setNodeValue(j, value); }
                } restore{n, j, v};

                n.setNodeValue(j, v + h);
                for (Size k = 0; k < rows.size(); ++k)
                    up[k] = value(rows[k]);

                n.setNodeValue(j, v - h);
                for (Size k = 0; k < rows.size(); ++k)
                    J[rows[k]][j - 1] =
                        (up[k] - value(rows[k])) / (2.0 * h);
            }
        }

        /*! Jacobian of continuous zero rates at node dates with respect to
            free stored nodes. Unsupported rows use numerical differences.
        */
        inline Matrix zeroNodeJacobian(
                const CurveJacobianNode& n,
                std::vector<bool>* analyticDerivatives = nullptr) {
            n.ensure();
            std::vector<Date> dates = n.nodeDates();
            Size rows = dates.size(), cols = n.numNodes();
            QL_REQUIRE(rows == cols,
                       "the curve has " << rows << " node dates but " <<
                       cols << " free nodes");

            Matrix J(rows, cols, 0.0);
            std::vector<bool> analytic(rows, false);
            std::vector<Real> row;
            auto zeroAt = [&](Size i) {
                Time t = n.curve->timeFromReference(dates[i]);
                DiscountFactor p = n.curve->discount(dates[i], true);
                QL_REQUIRE(t > 0.0 && p > 0.0,
                           "cannot calculate a continuous zero rate at node " << i);
                return -std::log(p) / t;
            };
            for (Size i = 0; i < rows; ++i) {
                Time t = n.curve->timeFromReference(dates[i]);
                DiscountFactor p = n.curve->discount(dates[i], true);
                QL_REQUIRE(t > 0.0 && p > 0.0,
                           "cannot calculate a continuous zero rate at node " << i);
                std::vector<std::pair<Date, Real>> sensitivity = {
                    {dates[i], -1.0 / (t * p)}};
                if (n.analyticRow(sensitivity, row)) {
                    QL_REQUIRE(row.size() == cols,
                               "analytical zero/node row has " << row.size() <<
                               " entries instead of " << cols);
                    std::copy(row.begin(), row.end(), J.row_begin(i));
                    analytic[i] = true;
                }
            }

            std::vector<Size> numericalRows;
            for (Size i = 0; i < rows; ++i)
                if (!analytic[i])
                    numericalRows.push_back(i);
            fillNumericalNodeRows(n, numericalRows, J, zeroAt);

            if (analyticDerivatives != nullptr)
                *analyticDerivatives = std::move(analytic);
            return J;
        }

        /*! Jacobian of free stored nodes with respect to continuous zero
            rates at node dates.
        */
        inline Matrix nodeZeroJacobian(
                const CurveJacobianNode& n,
                std::vector<bool>* analyticDerivatives = nullptr) {
            Matrix J = zeroNodeJacobian(n, analyticDerivatives);
            QL_REQUIRE(J.rows() == J.columns(),
                       "cannot invert a non-square zero/node Jacobian");
            return inverse(J);
        }

        //! implied-quote sensitivities of a curve's alive helpers, in row order
        /*! Hoist this out of loops over column curves: the result is
            identical for every column curve of the same row curve.
        */
    inline std::vector<QuoteSensitivities> aliveHelperSensitivities(const CurveJacobianNode& n) {
            auto helpers = n.aliveHelpers();
            std::vector<QuoteSensitivities> result;
            result.reserve(helpers.size());
            for (const auto& helper : helpers)
                result.push_back(helper->impliedQuoteSensitivitiesByCurve());
            return result;
        }

        /*! Jacobian of one curve's helper quotes with respect to another
            curve's nodes. The same-curve case uses the diagonal block.
            Unsupported rows are differentiated numerically. Analytical rows
            require complete dependency transforms. The optional
            rowSensitivities avoid recomputing helper sensitivities per
            column curve; entries must follow a's alive helpers.
        */
        inline Matrix curveCrossJacobian(const CurveJacobianNode& a,
                                         const CurveJacobianNode& b,
                                         const CurveCrossJacobianContext& context,
                                         std::vector<bool>* analyticEquations = nullptr,
                                         const std::vector<QuoteSensitivities>*
                                             rowSensitivities = nullptr) {
            a.ensure();
            b.ensure();
            auto helpers = a.aliveHelpers();
            Size rows = helpers.size();
            Size cols = b.numNodes();
            QL_REQUIRE(rowSensitivities == nullptr ||
                           rowSensitivities->size() == rows,
                       "cached helper sensitivities have " <<
                       (rowSensitivities == nullptr ? 0 : rowSensitivities->size()) <<
                       " rows instead of " << rows);
            Matrix J(rows, cols, 0.0);
            std::vector<bool> analytic(rows, false);

            std::vector<Real> row;
            for (Size i = 0; i < rows; ++i) {
                QuoteSensitivities fetched;
                if (rowSensitivities == nullptr)
                    fetched = helpers[i]->impliedQuoteSensitivitiesByCurve();
                const QuoteSensitivities& s = rowSensitivities != nullptr
                    ? (*rowSensitivities)[i] : fetched;
                if (!s.available || s.incomplete.count(b.id) != 0)
                    continue;
                bool allAccounted = true;
                for (const auto& [curve, entries] : s.sensitivities)
                    allAccounted = allAccounted && context.accountsFor(curve, b.id);
                for (const auto* curve : s.incomplete)
                    allAccounted = allAccounted && context.accountsFor(curve, b.id);
                if (!allAccounted)
                    continue;
                auto bucket = s.sensitivities.find(b.id);
                bool indirectDependence = false;
                bool transformed = true;
                std::vector<std::pair<Date, Real>> targetSensitivities;
                if (bucket != s.sensitivities.end())
                    targetSensitivities = bucket->second;
                for (const auto& [source, entries] : s.sensitivities) {
                    if (source == b.id ||
                        !context.dependsOn(source, b.id))
                        continue;
                    indirectDependence = true;
                    transformed = transformed &&
                                  context.propagate(
                                      source, b.id, entries, targetSensitivities);
                }
                for (const auto* source : s.incomplete)
                    if (source != b.id &&
                        context.dependsOn(source, b.id)) {
                        indirectDependence = true;
                        transformed = false;
                    }
                if (indirectDependence) {
                    if (transformed && b.analyticRow(targetSensitivities, row)) {
                        std::copy(row.begin(), row.end(), J.row_begin(i));
                        analytic[i] = true;
                    }
                } else if (bucket == s.sensitivities.end()) {
                    // no dependence on this curve
                    analytic[i] = true;
                } else if (b.analyticRow(bucket->second, row)) {
                    std::copy(row.begin(), row.end(), J.row_begin(i));
                    analytic[i] = true;
                }
            }

            std::vector<Size> numericalRows;
            for (Size i = 0; i < rows; ++i)
                if (!analytic[i])
                    numericalRows.push_back(i);
            fillNumericalNodeRows(
                b, numericalRows, J,
                [&](Size i) { return helpers[i]->impliedQuote(); });

            if (analyticEquations != nullptr)
                *analyticEquations = analytic;
            return J;
        }

        /*! Inverse Jacobian for a curve group.

              \f[ \sum_P J_{XP} \, dz_P = dq_X \f]

            Rows are nodes and columns are alive helpers in curve order.
            Optional offsets include the final total. A quote is analytical
            only when all its blocks are analytical.
        */
        inline Matrix curveGroupInverseJacobian(
                const std::vector<CurveJacobianNode>& curves,
                std::vector<Size>* rowOffsets = nullptr,
                std::vector<Size>* colOffsets = nullptr,
                std::vector<bool>* analyticEquations = nullptr,
                const CurveCrossJacobianContext& baseContext = {}) {
            Size n = curves.size();
            std::vector<Size> nodeOffset(n + 1, 0), quoteOffset(n + 1, 0);
            CurveCrossJacobianContext context = baseContext;
            for (Size i = 0; i < n; ++i) {
                curves[i].ensure();
                nodeOffset[i + 1] = nodeOffset[i] + curves[i].numNodes();
                quoteOffset[i + 1] = quoteOffset[i] + curves[i].aliveHelpers().size();
                context.addCurve(curves[i].id, curves[i].valueDependencies);
            }
            QL_REQUIRE(nodeOffset[n] == quoteOffset[n],
                       "cannot solve for the node/quote sensitivities: the "
                       "curves have " << nodeOffset[n] << " free nodes in "
                       "total but " << quoteOffset[n] << " alive helpers");

            Matrix J(quoteOffset[n], nodeOffset[n], 0.0);
            std::vector<bool> analytic(quoteOffset[n], true);
            for (Size i = 0; i < n; ++i) {
                std::vector<QuoteSensitivities> rowSensitivities =
                    aliveHelperSensitivities(curves[i]);
                for (Size j = 0; j < n; ++j) {
                    std::vector<bool> blockFlags;
                    Matrix block = curveCrossJacobian(curves[i], curves[j],
                                                      context, &blockFlags,
                                                      &rowSensitivities);
                    for (Size r = 0; r < block.rows(); ++r) {
                        std::copy(block.row_begin(r), block.row_end(r),
                                  J.row_begin(quoteOffset[i] + r) + nodeOffset[j]);
                        analytic[quoteOffset[i] + r] =
                            analytic[quoteOffset[i] + r] && blockFlags[r];
                    }
                }
            }

            if (rowOffsets != nullptr)
                *rowOffsets = nodeOffset;
            if (colOffsets != nullptr)
                *colOffsets = quoteOffset;
            if (analyticEquations != nullptr)
                *analyticEquations = analytic;
            return inverse(J);
        }

    }

}

#endif
