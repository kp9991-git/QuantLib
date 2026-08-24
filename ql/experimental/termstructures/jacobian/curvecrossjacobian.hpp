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
                n.aliveHelpers = [curve] {
                    curve->calculate();
                    Date firstDate = Traits::initialDate(curve.get());
                    std::vector<ext::shared_ptr<BootstrapHelper<YieldTermStructure>>> alive;
                    for (const auto& helper : curve->instruments_)
                        if (helper->pillarDate() > firstDate)
                            alive.push_back(helper);
                    return alive;
                };
                n.ownJacobian = [curve](std::vector<bool>* analyticRows) {
                    return curve->jacobian(analyticRows);
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
                    curve->jacobianCacheValid_ = false;
                    Traits::updateGuess(curve->data_, v, j);
                    curve->interpolation_.update();
                };
                return n;
            }
        };

        /*! Jacobian of one curve's helper quotes with respect to another
            curve's nodes. The same-curve case uses the diagonal block.
            Unsupported rows are differentiated numerically. Analytical rows
            require complete dependency transforms.
        */
        inline Matrix curveCrossJacobian(const CurveJacobianNode& a,
                                         const CurveJacobianNode& b,
                                         const CurveCrossJacobianContext& context,
                                         std::vector<bool>* analyticRows = nullptr) {
            a.ensure();
            b.ensure();
            auto helpers = a.aliveHelpers();
            Size rows = helpers.size();
            Size cols = b.numNodes();
            Matrix J(rows, cols, 0.0);
            std::vector<bool> analytic(rows, false);

            std::vector<Real> row;
            for (Size i = 0; i < rows; ++i) {
                auto s = helpers[i]->impliedQuoteSensitivitiesByCurve();
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
            if (!numericalRows.empty()) {
                std::vector<Real> up(numericalRows.size());
                for (Size j = 1; j <= cols; ++j) {
                    Real v = b.nodeValue(j);
                    Real h = 1.0e-6 * std::max(std::abs(v), 0.01);

                    // restore the node after success or failure
                    struct RestoreNode {  // NOLINT(cppcoreguidelines-special-member-functions)
                        const CurveJacobianNode& b;
                        Size j;
                        Real v;
                        ~RestoreNode() { b.setNodeValue(j, v); }
                    } restore{b, j, v};

                    b.setNodeValue(j, v + h);
                    for (Size k = 0; k < numericalRows.size(); ++k)
                        up[k] = helpers[numericalRows[k]]->impliedQuote();

                    b.setNodeValue(j, v - h);
                    for (Size k = 0; k < numericalRows.size(); ++k)
                        J[numericalRows[k]][j-1] =
                            (up[k] - helpers[numericalRows[k]]->impliedQuote())/(2.0*h);
                }
            }

            if (analyticRows != nullptr)
                *analyticRows = analytic;
            return J;
        }

        /*! Node sensitivities to all helper quotes for a curve group.

              \f[ \sum_P J_{XP} \, dz_P = dq_X \f]

            Rows are nodes and columns are alive helpers in curve order.
            Optional offsets include the final total. A quote is analytical
            only when all its blocks are analytical.
        */
        inline Matrix groupNodeQuoteJacobian(const std::vector<CurveJacobianNode>& curves,
                                             std::vector<Size>* rowOffsets = nullptr,
                                             std::vector<Size>* colOffsets = nullptr,
                                             std::vector<bool>* analyticRows = nullptr,
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
                for (Size j = 0; j < n; ++j) {
                    std::vector<bool> blockFlags;
                    Matrix block = curveCrossJacobian(curves[i], curves[j],
                                                      context, &blockFlags);
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
            if (analyticRows != nullptr)
                *analyticRows = analytic;
            return inverse(J);
        }

    }

}

#endif
