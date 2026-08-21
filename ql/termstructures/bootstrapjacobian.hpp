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

/*! \file bootstrapjacobian.hpp
    \brief shared implementation of the Jacobian of bootstrapped curves
*/

#ifndef quantlib_bootstrap_jacobian_hpp
#define quantlib_bootstrap_jacobian_hpp

#include <ql/experimental/termstructures/curvechainrulecalculator.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/matrix.hpp>
#include <ql/termstructures/bootstraphelper.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <algorithm>
#include <functional>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace QuantLib {

    namespace detail {

        // optional traits interface for analytical Jacobians
        template <class T, class Curve, class = void>
        constexpr bool hasSensitivityScale = false;

        template <class T, class Curve>
        constexpr bool hasSensitivityScale<
            T, Curve,
            std::void_t<decltype(T::sensitivityScale(
                Time(), std::declval<const Curve*>()))>> = true;

        template <class T, class = void>
        constexpr bool hasFirstDataPointFlag = false;

        template <class T>
        constexpr bool hasFirstDataPointFlag<
            T, std::void_t<decltype(T::firstDataPointTracksSecond)>> = true;

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

        //! type-erased bootstrapped curve used in cross-curve Jacobians
        /*! A null curve means the interface is unavailable. */
        struct CurveJacobianNode {
            ext::shared_ptr<YieldTermStructure> curve;
            CurveId id = nullptr;
            //! term structures used internally to turn nodes into curve values
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

        //! system curves and known dependent wrappers in a multi-curve group
        struct CurveJacobianGroup {
            std::vector<CurveJacobianNode> members;
            //! index of the curve that requested the group
            Size target = 0;
            std::set<const TermStructure*> dependents;
        };

        // access to curve internals needed for cross-curve Jacobians
        template <class Curve>
        struct BootstrapJacobianAccess;

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

        /*! Converts sensitivities to curve values \f$ v(t) \f$ into one
            row of node sensitivities. Returns false if traits, interpolation,
            or extrapolation are unsupported.
        */
        template <class Traits, class Curve>
        bool analyticNodeRow(const Curve* curve,
                             const std::vector<std::pair<Time, Real>>& sensitivities,
                             const std::vector<Time>& times,
                             const Interpolation& interpolation,
                             std::vector<Real>& row) {
            if constexpr (!hasSensitivityScale<Traits, Curve>) {
                return false;
            } else {
                bool firstTied = false;
                if constexpr (hasFirstDataPointFlag<Traits>)
                    firstTied = Traits::firstDataPointTracksSecond;

                Time tMax = times.back();
                row.resize(times.size() - 1);
                std::fill(row.begin(), row.end(), 0.0);
                for (const auto& [t, dQdP] : sensitivities) {
                    // curve extrapolation might differ from interpolation
                    if (t > tMax)
                        return false;
                    auto weights = interpolation.nodeWeights(t, true);
                    if (weights.empty())
                        return false;
                    Real scale = Traits::sensitivityScale(t, curve);
                    for (const auto& [j, w] : weights) {
                        if (j > 0)
                            row[j-1] += dQdP * scale * w;
                        else if (firstTied)
                            // data[0] tracks data[1]
                            row[0] += dQdP * scale * w;
                        // otherwise data[0] is fixed
                    }
                }
                return true;
            }
        }

        //! Jacobian of helper quotes with respect to curve nodes
        template <class Traits, class Curve>
        Matrix bootstrapJacobian(
                 const Curve* curve,
                 const std::vector<ext::shared_ptr<typename Traits::helper>>& instruments,
                 const std::vector<Time>& times,
                 std::vector<Real>& data,
                 Interpolation& interpolation,
                 bool curveHasJumps,
                 std::vector<bool>* analyticRows,
                 bool numericalFallback = true) {

            // alive helpers in curve order
            Date firstDate = Traits::initialDate(curve);
            std::vector<ext::shared_ptr<typename Traits::helper>> alive;
            for (const auto& helper : instruments)
                if (helper->pillarDate() > firstDate)
                    alive.push_back(helper);

            Size rows = alive.size();
            Size cols = times.size() - 1;
            Matrix J(rows, cols, 0.0);
            std::vector<bool> analytic(rows, false);

            if constexpr (hasSensitivityScale<Traits, Curve>) {
                // analytical weights do not include jumps
                if (!curveHasJumps) {
                    std::vector<Real> row(cols);
                    for (Size i = 0; i < rows; ++i) {
                        auto sensitivities = alive[i]->impliedQuoteSensitivities();
                        if (sensitivities.empty())
                            continue;
                        if (analyticNodeRow<Traits>(curve, sensitivities, times,
                                                    interpolation, row)) {
                            std::copy(row.begin(), row.end(), J.row_begin(i));
                            analytic[i] = true;
                        }
                    }
                }
            }

            // differentiate remaining rows one column at a time
            std::vector<Size> numericalRows;
            if (numericalFallback)
                for (Size i = 0; i < rows; ++i)
                    if (!analytic[i])
                        numericalRows.push_back(i);

            if (!numericalRows.empty()) {
                std::vector<Real> savedData = data;

                // restore element-wise to preserve interpolation iterators
                // the guard also restores after an exception
                struct RestoreData {  // NOLINT(cppcoreguidelines-special-member-functions)
                    std::vector<Real>& data;
                    const std::vector<Real>& saved;
                    Interpolation& interpolation;
                    ~RestoreData() {
                        std::copy(saved.begin(), saved.end(), data.begin());
                        interpolation.update();
                    }
                } restore{data, savedData, interpolation};

                std::vector<Real> up(numericalRows.size());
                for (Size j = 1; j <= cols; ++j) {
                    Real v = savedData[j];
                    Real h = 1.0e-6 * std::max(std::abs(v), 0.01);

                    Traits::updateGuess(data, v + h, j);
                    interpolation.update();
                    for (Size k = 0; k < numericalRows.size(); ++k)
                        up[k] = alive[numericalRows[k]]->impliedQuote();

                    Traits::updateGuess(data, v - h, j);
                    interpolation.update();
                    for (Size k = 0; k < numericalRows.size(); ++k)
                        J[numericalRows[k]][j-1] =
                            (up[k] - alive[numericalRows[k]]->impliedQuote())/(2.0*h);

                    std::copy(savedData.begin(), savedData.end(), data.begin());
                    interpolation.update();
                }
            }

            if (analyticRows != nullptr)
                *analyticRows = analytic;
            return J;
        }

        //! invert a square curve Jacobian
        inline Matrix inverseBootstrapJacobian(const Matrix& J) {
            QL_REQUIRE(J.rows() == J.columns(),
                       "cannot invert the Jacobian: the curve has " <<
                       J.columns() << " free nodes but only " << J.rows() <<
                       " alive helpers");
            return inverse(J);
        }

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
                    return analyticNodeRow<Traits>(curve.get(), sensitivities,
                                                   curve->times_, curve->interpolation_,
                                                   row);
                };
                n.nodeValue = [curve](Size j) {
                    curve->calculate();
                    return curve->data_[j];
                };
                n.setNodeValue = [curve](Size j, Real v) {
                    // preserve interpolation iterators
                    Traits::updateGuess(curve->data_, v, j);
                    curve->interpolation_.update();
                };
                return n;
            }
        };

        /*! Partial Jacobian of the first curve's helper quotes with respect
            to the second curve's nodes. Coincident curves use the own-curve
            Jacobian. Unsupported rows are differentiated numerically.

            Analytical rows require every referenced curve to be accounted
            for by the supplied context. Known dependent wrappers without an
            analytical dependency transform remain numerical.
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
                    // no direct dependence on this curve
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

        /*! Node sensitivities to all helper quotes, obtained by solving
            the differentiated bootstrap conditions

              \f[ \sum_P J_{XP} \, dz_P = dq_X \f]

            for the whole group. Rows are nodes and columns are alive helpers,
            both grouped in the given curve order. Optional offsets include a
            final total. A quote is analytical only if all its blocks are.

            Group members and their dependency edges are added to a copy of
            the supplied context.
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
