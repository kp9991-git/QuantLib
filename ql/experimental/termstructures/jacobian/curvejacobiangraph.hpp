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
#include <ql/experimental/termstructures/jacobian/curveriskpropagation.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <algorithm>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace QuantLib {

    //! cross-curve Jacobians for bootstrapped curves
    /*! Combines the bootstrap equations of registered curves. Rows follow
        alive helpers and columns follow free nodes in registration order.
    */
    class CurveJacobianGraph {
      public:
        explicit CurveJacobianGraph(bool errorOnIncomplete = false)
        : errorOnIncomplete_(errorOnIncomplete) {}

        /*! Whether every curve dependency reported by registered curves and
            their helpers is represented in the graph. Helpers without
            dependency metadata do not make the graph incomplete.
        */
        bool isComplete() const {
            for (const auto& node : nodes_) {
                if (!dependenciesAreRegistered(
                        node.valueDependencies.targets()))
                    return false;
                for (const auto& helper : node.aliveHelpers()) {
                    QuoteSensitivities sensitivities =
                        helper->impliedQuoteSensitivitiesByCurve();
                    for (const auto& sensitivity : sensitivities.sensitivities)
                        if (!isRegistered(sensitivity.first))
                            return false;
                    if (!dependenciesAreRegistered(sensitivities.incomplete))
                        return false;
                }
            }
            for (detail::CurveId derived : derivedIds_)
                if (!dependenciesAreRegistered(
                        derivedDependencies_.targets(derived)))
                    return false;
            return true;
        }

        /*! Register a calibrated curve or a supported derived curve.
            Calibrated curves contribute blocks. Derived curves only record
            dependencies.
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

        /*! Partial Jacobian of the first curve's helper quotes with respect
            to the second curve's nodes. Other curve nodes are fixed.
        */
        Matrix crossJacobian(const YieldTermStructure& of,
                             const YieldTermStructure& withRespectTo,
            std::vector<bool>* analyticEquations = nullptr) const {
            requireComplete();
            detail::CurveCrossJacobianContext context = jacobianContext();
            return detail::curveCrossJacobian(node(of), node(withRespectTo),
                                              context, analyticEquations);
        }

        /*! Block of the inverse Jacobian for the registered curve system.
            Rows are the first curve's nodes and columns are the second
            curve's helper quotes. All registered dependencies are included.
            The optional flags describe the second curve's helper equations.
        */
        Matrix inverseJacobian(const YieldTermStructure& of,
                               const YieldTermStructure& withRespectTo,
                               std::vector<bool>* analyticEquations = nullptr) const {
            requireComplete();
            Size a = index(of), b = index(withRespectTo);
            detail::CurveCrossJacobianContext context = jacobianContext(false);
            detail::CurveJacobianBlocks blocks =
                detail::curveJacobianBlocks(nodes_, context);
            Matrix block = detail::inverseCurveJacobianBlock(blocks, a, b);
            if (analyticEquations != nullptr)
                *analyticEquations = blocks.analyticQuotes[b];
            return block;
        }

        /*! Block of the dense inverse Jacobian for the registered curve
            system. This is the reference implementation of inverseJacobian().
        */
        Matrix inverseJacobianDense(
                const YieldTermStructure& of,
                const YieldTermStructure& withRespectTo,
                std::vector<bool>* analyticEquations = nullptr) const {
            requireComplete();
            Size a = index(of), b = index(withRespectTo);
            detail::CurveGroupInverse dense =
                detail::curveGroupInverseJacobian(nodes_, jacobianContext(false));
            const std::vector<Size>& rows = dense.nodeOffset;
            const std::vector<Size>& cols = dense.quoteOffset;

            Matrix block(rows[a + 1] - rows[a], cols[b + 1] - cols[b]);
            for (Size i = 0; i < block.rows(); ++i)
                std::copy(
                    dense.inverse.row_begin(rows[a] + i) + cols[b],
                    dense.inverse.row_begin(rows[a] + i) + cols[b + 1],
                    block.row_begin(i));

            if (analyticEquations != nullptr)
                *analyticEquations = std::vector<bool>(
                    dense.analyticEquations.begin() + cols[b],
                    dense.analyticEquations.begin() + cols[b + 1]);
            return block;
        }

        /*! Jacobian of continuously compounded zero rates at a curve's node
            dates with respect to its free stored nodes.
        */
        Matrix zeroNodeJacobian(
                const YieldTermStructure& curve,
                std::vector<bool>* analyticDerivatives = nullptr) const {
            requireComplete();
            return detail::zeroNodeJacobian(node(curve), analyticDerivatives);
        }

        /*! Jacobian of a curve's free stored nodes with respect to
            continuously compounded zero rates at its node dates.
        */
        Matrix nodeZeroJacobian(
                const YieldTermStructure& curve,
                std::vector<bool>* analyticDerivatives = nullptr) const {
            requireComplete();
            return detail::nodeZeroJacobian(node(curve), analyticDerivatives);
        }

        /*! Convert node risk to par risk using the dense inverse.
            This is the reference implementation of parRisk().
        */
        std::map<const YieldTermStructure*, Array> parRiskDense(
            const std::map<const YieldTermStructure*, Array>& nodeRisk,
                     std::vector<bool>* analyticEquations = nullptr) const {
            requireComplete();
            detail::CurveGroupInverse dense =
                detail::curveGroupInverseJacobian(nodes_, jacobianContext(false));
            const std::vector<Size>& rows = dense.nodeOffset;
            const std::vector<Size>& cols = dense.quoteOffset;
            if (analyticEquations != nullptr)
                *analyticEquations = std::move(dense.analyticEquations);

            // s is the stacked node risk and r = transpose(S) * s
            Array s(dense.inverse.rows(), 0.0);
            for (const auto& [curve, risk] : nodeRisk) {
                Size a = index(*curve);
                QL_REQUIRE(risk.size() == rows[a+1] - rows[a],
                           "node risk size (" << risk.size() <<
                           ") does not match the number of curve nodes (" <<
                           rows[a+1] - rows[a] << ")");
                std::copy(risk.begin(), risk.end(), s.begin() + rows[a]);
            }
            Array r = transpose(dense.inverse)*s;

            std::map<const YieldTermStructure*, Array> result;
            for (Size b = 0; b < nodes_.size(); ++b)
                result[nodes_[b].curve.get()] =
                    Array(r.begin() + cols[b], r.begin() + cols[b+1]);
            return result;
        }

        /*! Propagate direct node risk over the dependency graph.
            Par risk solves \f$ J^T r=s \f$ in each dependency component.
            Zero risk fixes each curve while the others rebootstrap. In an
            acyclic component it satisfies
            \f$ A_{bb}^T r_b=z_b \f$. A cyclic component requires a reduced
            solve and need not satisfy this identity. Either output may be
            null. Input arrays must match node counts.
        */
        void propagateNodeRisk(
                const std::map<const YieldTermStructure*, Array>& nodeRisk,
                std::map<const YieldTermStructure*, Array>* zeroRisk,
                std::map<const YieldTermStructure*, Array>* parRisk,
                std::vector<bool>* analyticEquations = nullptr) const {
            requireComplete();
            detail::CurveCrossJacobianContext context = jacobianContext(false);
            detail::CurveJacobianBlocks blocks =
                detail::curveJacobianBlocks(nodes_, context);

            std::vector<Array> direct(nodes_.size());
            for (const auto& [curve, risk] : nodeRisk) {
                Size a = index(*curve);
                QL_REQUIRE(risk.size() == blocks.numNodes(a),
                           "node risk size (" << risk.size() <<
                           ") does not match the number of curve nodes (" <<
                           blocks.numNodes(a) << ")");
                direct[a] = risk;
            }

            detail::CurveRiskPropagation propagated =
                detail::propagateCurveNodeRisk(blocks, direct,
                                               zeroRisk != nullptr);

            for (Size b = 0; b < nodes_.size(); ++b) {
                if (zeroRisk != nullptr)
                    (*zeroRisk)[nodes_[b].curve.get()] = propagated.nodeRisk[b];
                if (parRisk != nullptr)
                    (*parRisk)[nodes_[b].curve.get()] = propagated.quoteRisk[b];
            }

            if (analyticEquations != nullptr) {
                std::vector<bool> flat;
                flat.reserve(blocks.quoteOffset.back());
                for (const auto& curveFlags : blocks.analyticQuotes)
                    flat.insert(flat.end(),
                                curveFlags.begin(), curveFlags.end());
                *analyticEquations = std::move(flat);
            }
        }

        /*! Convert direct node risk to node risk for all registered curves,
            following the dependency graph. Input keys must be registered and
            arrays must match node counts.
        */
        std::map<const YieldTermStructure*, Array> zeroRisk(
            const std::map<const YieldTermStructure*, Array>& nodeRisk,
                 std::vector<bool>* analyticEquations = nullptr) const {
            std::map<const YieldTermStructure*, Array> result;
            propagateNodeRisk(nodeRisk, &result, nullptr, analyticEquations);
            return result;
        }

        /*! Convert node risk to par-instrument risk for all registered curves.
            Input keys must be registered and arrays must match node counts.
        */
        std::map<const YieldTermStructure*, Array> parRisk(
            const std::map<const YieldTermStructure*, Array>& nodeRisk,
                std::vector<bool>* analyticEquations = nullptr) const {
            std::map<const YieldTermStructure*, Array> result;
            propagateNodeRisk(nodeRisk, nullptr, &result, analyticEquations);
            return result;
        }

      private:
        void requireComplete() const {
            QL_REQUIRE(!errorOnIncomplete_ || isComplete(),
                       "incomplete curve Jacobian graph");
        }

        //! whether the curve was added, as a calibrated or derived curve
        bool isRegistered(detail::CurveId id) const {
            if (derivedIds_.count(id) != 0)
                return true;
            return std::any_of(
                nodes_.begin(), nodes_.end(),
                [id](const detail::CurveJacobianNode& node) {
                    return node.id == id;
                });
        }

        template <class Container>
        bool dependenciesAreRegistered(const Container& dependencies) const {
            return std::all_of(
                dependencies.begin(), dependencies.end(),
                [this](detail::CurveId id) { return isRegistered(id); });
        }

        void addDerivedCurve(
                const ext::shared_ptr<YieldTermStructure>& curve,
                const ext::shared_ptr<YieldTermStructure>& dependency) {
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

        detail::CurveCrossJacobianContext jacobianContext(bool includeAccountedCurves = true) const {
            detail::CurveCrossJacobianContext result;
            result.assumeUnlistedCurvesIndependent();
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
        bool errorOnIncomplete_;
    };

}

#endif
