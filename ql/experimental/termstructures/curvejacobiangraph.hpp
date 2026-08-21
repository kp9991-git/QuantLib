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

/*! \file curvejacobiangraph.hpp
    \brief cross-curve Jacobians of a set of bootstrapped curves
*/

#ifndef quantlib_experimental_curve_jacobian_graph_hpp
#define quantlib_experimental_curve_jacobian_graph_hpp

#include <ql/math/array.hpp>
#include <ql/math/matrix.hpp>
#include <ql/termstructures/bootstrapjacobian.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace QuantLib {

    //! cross-curve Jacobians for bootstrapped curves
    /*! Combines the bootstrap equations of registered curves, including
        acyclic dependencies and curves bootstrapped jointly by MultiCurve.
        Rows follow alive helpers and columns follow free curve nodes in
        registration order.
    */
    class CurveJacobianGraph {
      public:
        /*! Register a calibrated curve or inspect a supported derived curve.

            Calibrated curves contribute nodes and helper quotes to the
            Jacobian system.  Non-calibrated curves exposing baseCurve() are
            recorded as dependency bridges and do not contribute a block.
        */
        template <class Curve>
        void add(const ext::shared_ptr<Curve>& curve) {
            QL_REQUIRE(curve, "null curve");
            if constexpr (detail::supportsCurveJacobianNode<Curve>) {
                detail::CurveJacobianNode n =
                    detail::BootstrapJacobianAccess<Curve>::makeNode(curve);
                for (auto& existing : nodes_) {
                    if (existing.id == n.id) {
                        existing = std::move(n);
                        return;
                    }
                }
                nodes_.push_back(std::move(n));
            } else if constexpr (detail::hasBaseCurveHandle<Curve>) {
                QL_REQUIRE(!curve->baseCurve().empty(),
                           "derived curve has an empty base-curve handle");
                addDerivedCurve(curve, curve->baseCurve().currentLink());
            } else {
                QL_FAIL("curve type is not supported by CurveJacobianGraph::add; "
                        "supported curves must provide bootstrap Jacobians or "
                        "expose baseCurve()");
            }
        }

        /*! Validate all dependencies reported by the registered helpers.

            If requireAnalyticMetadata is true, helpers that do not expose
            QuoteSensitivities are rejected as well.  Otherwise they remain
            eligible for numerical differentiation, but dependencies reported
            by helpers that do expose metadata are always checked.
        */
        void validateDependencies(bool requireAnalyticMetadata = false) const {
            std::set<const TermStructure*> declared = accountedCurves();
            declared.insert(derivedIds_.begin(), derivedIds_.end());
            detail::CurveChainRuleCalculator chainRule = chainRuleCalculator();

            const std::set<const TermStructure*> terminal = accountedCurves();
            std::function<void(const TermStructure*, std::set<const TermStructure*>&)>
                validateDependentPath;
            validateDependentPath = [&](const TermStructure* dependency,
                                        std::set<const TermStructure*>& visiting) {
                if (terminal.count(dependency) != 0)
                    return;
                std::vector<detail::CurveId> outgoing =
                    chainRule.targets(dependency);
                QL_REQUIRE(!outgoing.empty(),
                           "a derived curve dependency was not added to "
                           "the curve Jacobian graph");
                QL_REQUIRE(visiting.insert(dependency).second,
                           "cyclic dependency between unregistered curve wrappers");
                for (const auto* target : outgoing)
                    validateDependentPath(target, visiting);
                visiting.erase(dependency);
            };
            for (const auto& wrapper : derivedCurves_) {
                const auto* wrapperId =
                    static_cast<const TermStructure*>(wrapper.get());
                std::set<const TermStructure*> visiting;
                validateDependentPath(wrapperId, visiting);
                for (const auto* target : chainRule.targets(wrapperId))
                    QL_REQUIRE(declared.count(target) != 0,
                               "a derived curve references an undeclared curve");
            }

            for (Size c = 0; c < nodes_.size(); ++c) {
                nodes_[c].ensure();
                for (const auto* target : nodes_[c].valueDependencies.targets()) {
                    QL_REQUIRE(target != nullptr &&
                                   declared.count(target) != 0,
                               "registered curve " << c <<
                               " uses a base curve that was not added or "
                               "discovered through a supported derived curve");
                    std::set<const TermStructure*> visiting;
                    validateDependentPath(target, visiting);
                }
                auto helpers = nodes_[c].aliveHelpers();
                for (Size h = 0; h < helpers.size(); ++h) {
                    QuoteSensitivities s =
                        helpers[h]->impliedQuoteSensitivitiesByCurve();
                    if (!s.available) {
                        QL_REQUIRE(!requireAnalyticMetadata,
                                   "helper " << h << " of registered curve " << c <<
                                   " does not expose dependency metadata");
                        continue;
                    }
                    for (const auto& [dependency, entries] : s.sensitivities) {
                        QL_REQUIRE(dependency != nullptr &&
                                       declared.count(dependency) != 0,
                                   "helper " << h << " of registered curve " << c <<
                                   " references a curve that was not added or "
                                   "discovered through a supported derived curve");
                    }
                    for (const auto* dependency : s.incomplete) {
                        QL_REQUIRE(dependency != nullptr &&
                                       declared.count(dependency) != 0,
                                   "helper " << h << " of registered curve " << c <<
                                   " has an incomplete contribution from a curve "
                                   "that was not added or discovered through a "
                                   "supported derived curve");
                    }
                }
            }
        }

        /*! Partial Jacobian of the first curve's helper quotes with respect
            to the second curve's nodes. Other curve nodes are fixed.
        */
        Matrix crossJacobian(const YieldTermStructure& of,
                             const YieldTermStructure& withRespectTo,
            std::vector<bool>* analyticRows = nullptr) const {
            validateDependencies();
            detail::CurveCrossJacobianContext context = jacobianContext();
            return detail::curveCrossJacobian(node(of), node(withRespectTo),
                                              context, analyticRows);
        }

        /*! Jacobian of the first curve's nodes with respect to the second
            curve's helper quotes, including all registered dependencies.
        */
        Matrix nodeQuoteJacobian(const YieldTermStructure& of,
                                 const YieldTermStructure& withRespectTo,
                                 std::vector<bool>* analyticRows = nullptr) const {
            validateDependencies();
            Size a = index(of), b = index(withRespectTo);
            std::vector<Size> rowOffsets, colOffsets;
            std::vector<bool> allAnalytic;
            detail::CurveCrossJacobianContext context = jacobianContext(false);
            Matrix S = detail::groupNodeQuoteJacobian(nodes_, &rowOffsets, &colOffsets,
                                                      &allAnalytic, context);
            Matrix block(rowOffsets[a+1] - rowOffsets[a],
                         colOffsets[b+1] - colOffsets[b]);
            for (Size i = 0; i < block.rows(); ++i)
                std::copy(S.row_begin(rowOffsets[a] + i) + colOffsets[b],
                          S.row_begin(rowOffsets[a] + i) + colOffsets[b+1],
                          block.row_begin(i));
            if (analyticRows != nullptr)
                *analyticRows = std::vector<bool>(allAnalytic.begin() + colOffsets[b],
                                                  allAnalytic.begin() + colOffsets[b+1]);
            return block;
        }

        /*! Convert node risk to par-instrument risk for all registered curves.
            Input keys must be registered and arrays must match node counts.
        */
        std::map<const YieldTermStructure*, Array>
        parRisk(const std::map<const YieldTermStructure*, Array>& nodeRisk,
                std::vector<bool>* analyticRows = nullptr) const {
            validateDependencies();
            std::vector<Size> rowOffsets, colOffsets;
            std::vector<bool> allAnalytic;
            detail::CurveCrossJacobianContext context = jacobianContext(false);
            Matrix S = detail::groupNodeQuoteJacobian(nodes_, &rowOffsets, &colOffsets,
                                                      &allAnalytic, context);
            if (analyticRows != nullptr)
                *analyticRows = std::move(allAnalytic);

            // r = transpose(S) * s, with s the stacked node risk
            Array s(S.rows(), 0.0);
            for (const auto& [curve, risk] : nodeRisk) {
                Size a = index(*curve);
                QL_REQUIRE(risk.size() == rowOffsets[a+1] - rowOffsets[a],
                           "node risk size (" << risk.size() <<
                           ") does not match the number of curve nodes (" <<
                           rowOffsets[a+1] - rowOffsets[a] << ")");
                std::copy(risk.begin(), risk.end(), s.begin() + rowOffsets[a]);
            }
            Array r = transpose(S)*s;

            std::map<const YieldTermStructure*, Array> result;
            for (Size b = 0; b < nodes_.size(); ++b)
                result[nodes_[b].curve.get()] =
                    Array(r.begin() + colOffsets[b], r.begin() + colOffsets[b+1]);
            return result;
        }

      private:
        void addDerivedCurve(
                const ext::shared_ptr<const YieldTermStructure>& curve,
                const ext::shared_ptr<const YieldTermStructure>& dependency) {
            QL_REQUIRE(dependency, "null base curve of derived curve");
            QL_REQUIRE(dependency.get() != curve.get(),
                       "a derived curve cannot depend on itself");
            const auto* id = static_cast<const TermStructure*>(curve.get());
            derivedIds_.insert(id);
            derivedDependencies_.add(
                id, static_cast<const TermStructure*>(dependency.get()));
            for (const auto& existing : derivedCurves_)
                if (existing.get() == curve.get())
                    return;
            derivedCurves_.push_back(curve);
        }

        Size index(const YieldTermStructure& curve) const {
            for (Size i = 0; i < nodes_.size(); ++i)
                if (nodes_[i].curve.get() == &curve)
                    return i;
            QL_FAIL("the given curve was not added to the graph");
        }

        const detail::CurveJacobianNode& node(const YieldTermStructure& curve) const {
            return nodes_[index(curve)];
        }

        std::set<const TermStructure*> accountedCurves() const {
            std::set<const TermStructure*> accounted;
            for (const auto& n : nodes_)
                accounted.insert(n.id);
            return accounted;
        }

        detail::CurveChainRuleCalculator chainRuleCalculator() const {
            auto result = derivedDependencies_;
            for (const auto& node : nodes_)
                result.add(node.id, node.valueDependencies);
            return result;
        }

        detail::CurveCrossJacobianContext
        jacobianContext(bool includeAccountedCurves = true) const {
            detail::CurveCrossJacobianContext result;
            result.addDependencies(derivedDependencies_);
            result.addNumericallyPropagatedCurves(derivedIds_);
            for (const auto& node : nodes_) {
                if (includeAccountedCurves)
                    result.addCurve(node.id, node.valueDependencies);
                else
                    result.addDependencies(node.id, node.valueDependencies);
            }
            return result;
        }

        std::vector<detail::CurveJacobianNode> nodes_;
        std::vector<ext::shared_ptr<const YieldTermStructure>> derivedCurves_;
        std::set<const TermStructure*> derivedIds_;
        detail::CurveChainRuleCalculator derivedDependencies_;
    };

}

#endif
