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

/*! \file curvejacobiangraph.hpp
    \brief cross-curve Jacobians of a set of bootstrapped curves
*/

#ifndef quantlib_curve_jacobian_graph_hpp
#define quantlib_curve_jacobian_graph_hpp

#include <ql/math/array.hpp>
#include <ql/math/matrix.hpp>
#include <ql/termstructures/bootstrapjacobian.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
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
        //! register or replace a curve
        template <class Curve>
        void add(const ext::shared_ptr<Curve>& curve) {
            detail::CurveJacobianNode n =
                detail::BootstrapJacobianAccess<Curve>::makeNode(curve);
            for (auto& existing : nodes_) {
                if (existing.id == n.id) {
                    existing = std::move(n);
                    return;
                }
            }
            nodes_.push_back(std::move(n));
        }

        /*! Mark an exogenous curve as independent of registered curves.
            Undeclared references use numerical differentiation to preserve
            dependencies through wrappers such as spreaded curves.
        */
        void addConstant(const ext::shared_ptr<const YieldTermStructure>& curve) {
            constants_.push_back(curve);
            constantIds_.insert(static_cast<const TermStructure*>(curve.get()));
        }

        /*! Partial Jacobian of the first curve's helper quotes with respect
            to the second curve's nodes. Other curve nodes are fixed.
        */
        Matrix crossJacobian(const YieldTermStructure& of,
                             const YieldTermStructure& withRespectTo,
                             std::vector<bool>* analyticRows = nullptr) const {
            std::set<const TermStructure*> accounted = accountedCurves();
            return detail::curveCrossJacobian(node(of), node(withRespectTo),
                                              analyticRows, &accounted);
        }

        /*! Jacobian of the first curve's nodes with respect to the second
            curve's helper quotes, including all registered dependencies.
        */
        Matrix nodeQuoteJacobian(const YieldTermStructure& of,
                                 const YieldTermStructure& withRespectTo) const {
            Size a = index(of), b = index(withRespectTo);
            std::vector<Size> rowOffsets, colOffsets;
            Matrix S = detail::groupNodeQuoteJacobian(nodes_, &rowOffsets, &colOffsets,
                                                      nullptr, &constantIds_);
            Matrix block(rowOffsets[a+1] - rowOffsets[a],
                         colOffsets[b+1] - colOffsets[b]);
            for (Size i = 0; i < block.rows(); ++i)
                std::copy(S.row_begin(rowOffsets[a] + i) + colOffsets[b],
                          S.row_begin(rowOffsets[a] + i) + colOffsets[b+1],
                          block.row_begin(i));
            return block;
        }

        /*! Convert node risk to helper-quote risk for all registered curves.
            Input keys must be registered and arrays must match node counts.
        */
        std::map<const YieldTermStructure*, Array>
        quoteRisk(const std::map<const YieldTermStructure*, Array>& nodeRisk) const {
            std::vector<Size> rowOffsets, colOffsets;
            Matrix S = detail::groupNodeQuoteJacobian(nodes_, &rowOffsets, &colOffsets,
                                                      nullptr, &constantIds_);

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
            std::set<const TermStructure*> accounted = constantIds_;
            for (const auto& n : nodes_)
                accounted.insert(n.id);
            return accounted;
        }

        std::vector<detail::CurveJacobianNode> nodes_;
        std::vector<ext::shared_ptr<const YieldTermStructure>> constants_;
        std::set<const TermStructure*> constantIds_;
    };

}

#endif
