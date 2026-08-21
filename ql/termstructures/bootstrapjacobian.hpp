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

/*! \file bootstrapjacobian.hpp
    \brief shared implementation of the Jacobian of bootstrapped curves
*/

#ifndef quantlib_bootstrap_jacobian_hpp
#define quantlib_bootstrap_jacobian_hpp

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

        // detection of the optional traits interface used for
        // analytical Jacobians
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

        //! type-erased view of a bootstrapped curve used to assemble
        //! cross-curve Jacobians; a null curve means not available
        struct CurveJacobianNode {
            ext::shared_ptr<YieldTermStructure> curve;
            const TermStructure* id = nullptr;
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

        //! the members of a multi-curve group, and the curves known to
        //! be functions of them without being system variables (e.g.
        //! spreaded curves over a member)
        struct CurveJacobianGroup {
            std::vector<CurveJacobianNode> members;
            std::set<const TermStructure*> dependents;
        };

        // grants CurveJacobianGraph and GlobalBootstrap access to the
        // curve internals needed to assemble cross-curve Jacobians
        template <class Curve>
        struct BootstrapJacobianAccess;

        // detection of bootstrap classes able to report the members of
        // a multi-curve group
        template <class B, class = void>
        constexpr bool hasJacobianGroup = false;

        template <class B>
        constexpr bool hasJacobianGroup<
            B, std::void_t<decltype(std::declval<const B&>().jacobianGroup())>> = true;

        // detection of curves supported by BootstrapJacobianAccess:
        // bootstrapped yield curves providing the jacobian() interface
        template <class Curve, class = void>
        constexpr bool supportsCurveJacobianNode = false;

        template <class Curve>
        constexpr bool supportsCurveJacobianNode<
            Curve, std::void_t<decltype(std::declval<const Curve&>().jacobian(
                std::declval<std::vector<bool>*>()))>> =
            std::is_base_of_v<YieldTermStructure, Curve>;

        /*! Converts sensitivities of a quantity to the curve values
            \f$ v(t) \f$ (discount factors, for yield curves) into
            sensitivities to the curve node values, filling one row of a
            Jacobian.  Returns false when the analytical machinery does
            not cover the curve (no traits support, node weights not
            implemented for the interpolation, or times beyond the last
            node, where the curve might extrapolate differently than the
            interpolation does).
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
                    // beyond the last node the curve might extrapolate
                    // differently than the interpolation does
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
                            // data[0] is kept equal to data[1]
                            row[0] += dQdP * scale * w;
                        // otherwise data[0] is fixed and
                        // contributes nothing
                    }
                }
                return true;
            }
        }

        /*! Jacobian of the implied quotes of the given helpers with
            respect to the curve node values.
        */
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

            // alive helpers, in the order stored in the curve
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
                // jumps are applied on top of the interpolated values,
                // and the analytical machinery does not account for them
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

            // numerical fallback for the remaining rows, one column at a time
            std::vector<Size> numericalRows;
            if (numericalFallback)
                for (Size i = 0; i < rows; ++i)
                    if (!analytic[i])
                        numericalRows.push_back(i);

            if (!numericalRows.empty()) {
                std::vector<Real> savedData = data;

                // restore the data even when an evaluation throws;
                // element-wise, because assigning the vector would
                // invalidate the iterators stored in the interpolation
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

        //! inverts a curve Jacobian, checking that it is square
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
                    // jumps are applied on top of the interpolated values,
                    // and the analytical machinery does not account for them
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
                    // single-element writes: the iterators stored in the
                    // interpolation stay valid
                    Traits::updateGuess(curve->data_, v, j);
                    curve->interpolation_.update();
                };
                return n;
            }
        };

        /*! Jacobian of the implied quotes of the first curve's helpers
            with respect to the nodes of the second curve, all other
            curves' nodes being kept fixed; the curves' own Jacobian
            when they coincide.  Rows that cannot be computed
            analytically fall back to numerical differentiation, bumping
            the second curve's nodes directly so that no bootstrap is
            triggered and the result remains a partial derivative.

            A row can be computed analytically only when every curve the
            helper references is accounted for: the curve being bumped,
            one of the curves in fixedCurves (other system variables, or
            curves known to be constant), or --- when
            unknownCurvesAreFixed is true --- any curve not listed in
            dependentCurves.  Rows referencing unaccounted curves fall
            back to numerical differentiation, which remains correct
            when those curves are functions of the bumped curve (e.g.
            spreaded or implied term structures computing off it on
            demand), while the analytical shortcut would wrongly treat
            them as constants.
        */
        inline Matrix curveCrossJacobian(const CurveJacobianNode& a,
                                         const CurveJacobianNode& b,
                                         std::vector<bool>* analyticRows = nullptr,
                                         const std::set<const TermStructure*>* fixedCurves = nullptr,
                                         const std::set<const TermStructure*>* dependentCurves = nullptr,
                                         bool unknownCurvesAreFixed = false) {
            if (a.id == b.id)
                return a.ownJacobian(analyticRows);

            a.ensure();
            b.ensure();
            auto helpers = a.aliveHelpers();
            Size rows = helpers.size();
            Size cols = b.numNodes();
            Matrix J(rows, cols, 0.0);
            std::vector<bool> analytic(rows, false);

            auto accounted = [&](const TermStructure* curve) {
                if (curve == b.id)
                    return true;
                if (fixedCurves != nullptr && fixedCurves->count(curve) != 0)
                    return true;
                return unknownCurvesAreFixed &&
                       (dependentCurves == nullptr || dependentCurves->count(curve) == 0);
            };

            std::vector<Real> row;
            for (Size i = 0; i < rows; ++i) {
                auto s = helpers[i]->impliedQuoteSensitivitiesByCurve();
                if (!s.available || s.incomplete.count(b.id) != 0)
                    continue;
                bool allAccounted = true;
                for (const auto& [curve, entries] : s.sensitivities)
                    allAccounted = allAccounted && accounted(curve);
                for (const auto* curve : s.incomplete)
                    allAccounted = allAccounted && accounted(curve);
                if (!allAccounted)
                    continue;
                auto bucket = s.sensitivities.find(b.id);
                if (bucket == s.sensitivities.end()) {
                    // the helper does not depend on the curve
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

                    // restore the node even when an evaluation throws
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

        /*! Sensitivities of the node values of all given curves with
            respect to the quotes of all their helpers, obtained by
            solving the differentiated bootstrap conditions

              \f[ \sum_P J_{XP} \, dz_P = dq_X \f]

            of the whole group at once; co-dependent (cyclical) groups
            are handled naturally.  Rows follow the curve nodes and
            columns the alive helpers, both grouped by curve in the
            given order; the offsets of each curve's first row and
            column are returned in the optional arguments, with a final
            entry equal to the total size.  The optional flag vector
            reports, for each quote of the group, whether all its blocks
            were computed analytically.

            The given curves are the system variables and are accounted
            for automatically; the remaining arguments are forwarded to
            curveCrossJacobian() and control how curves outside the
            system are treated.
        */
        inline Matrix groupNodeQuoteJacobian(const std::vector<CurveJacobianNode>& curves,
                                             std::vector<Size>* rowOffsets = nullptr,
                                             std::vector<Size>* colOffsets = nullptr,
                                             std::vector<bool>* analyticRows = nullptr,
                                             const std::set<const TermStructure*>* fixedCurves = nullptr,
                                             const std::set<const TermStructure*>* dependentCurves = nullptr,
                                             bool unknownCurvesAreFixed = false) {
            Size n = curves.size();
            std::vector<Size> nodeOffset(n + 1, 0), quoteOffset(n + 1, 0);
            std::set<const TermStructure*> fixedIds;
            if (fixedCurves != nullptr)
                fixedIds = *fixedCurves;
            for (Size i = 0; i < n; ++i) {
                curves[i].ensure();
                nodeOffset[i + 1] = nodeOffset[i] + curves[i].numNodes();
                quoteOffset[i + 1] = quoteOffset[i] + curves[i].aliveHelpers().size();
                fixedIds.insert(curves[i].id);
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
                                                      &blockFlags, &fixedIds,
                                                      dependentCurves,
                                                      unknownCurvesAreFixed);
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
