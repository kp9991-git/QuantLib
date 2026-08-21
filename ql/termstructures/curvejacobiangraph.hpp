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

    //! cross-curve Jacobians of a set of bootstrapped curves
    /*! Given a set of bootstrapped curves depending on each other (for
        instance, a projection curve built with an exogenous discount
        curve, or a cross-currency curve built with exogenous collateral
        and forecast curves, or co-dependent curves bootstrapped jointly
        through MultiCurve), this class provides the Jacobians of the
        helper quotes and curve nodes across the whole set, so that a
        sensitivity with respect to the nodes of any curve can be
        propagated back through the curve dependencies and expressed as
        a sensitivity with respect to the quoted instruments of any
        curve.

        Rows and columns follow the order of the alive helpers and of
        the curve nodes (excluding the value at the reference date,
        which is not a free variable) as stored in each curve.

        \par Example
        \code
        CurveJacobianGraph graph;
        graph.add(discountCurve);   // e.g. OIS
        graph.add(projectionCurve); // e.g. built with exogenous OIS discounting

        // node risk of a trade, e.g. from bumping the curve nodes
        std::map<const YieldTermStructure*, Array> nodeRisk = {
            {projectionCurve.get(), dNPVdNodes}};
        // quote risk on all registered curves
        auto risk = graph.quoteRisk(nodeRisk);
        \endcode
    */
    class CurveJacobianGraph {
      public:
        //! registers a curve; adding the same curve again replaces its entry
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

        /*! Declares a curve that the helpers may reference to be
            constant, i.e., independent of every registered curve.

            Helpers referencing curves that are neither registered nor
            declared constant have their rows computed by numerical
            differentiation instead of the analytical formulas: such
            curves might be functions of a registered curve (e.g. a
            spreaded or implied term structure over it), which the
            analytical shortcut would wrongly treat as fixed, while the
            numerical fallback prices through them correctly.  Declaring
            genuinely exogenous curves recovers the analytical rows.
        */
        void addConstant(const ext::shared_ptr<const YieldTermStructure>& curve) {
            constants_.push_back(curve);
            constantIds_.insert(static_cast<const TermStructure*>(curve.get()));
        }

        /*! Jacobian of the implied quotes of the first curve's helpers
            with respect to the nodes of the second curve, all other
            curves' nodes being kept fixed.  When the two curves
            coincide, this is the curve's own Jacobian.  Rows that could
            not be computed analytically fall back to numerical
            differentiation; the flag vector, if given, reports which
            rows were analytical.
        */
        Matrix crossJacobian(const YieldTermStructure& of,
                             const YieldTermStructure& withRespectTo,
                             std::vector<bool>* analyticRows = nullptr) const {
            std::set<const TermStructure*> accounted = accountedCurves();
            return detail::curveCrossJacobian(node(of), node(withRespectTo),
                                              analyticRows, &accounted);
        }

        /*! Jacobian of the nodes of the first curve with respect to the
            quotes of the second curve's helpers, obtained by solving
            the differentiated bootstrap conditions of all registered
            curves at once.  Co-dependent (cyclical) sets of curves are
            handled naturally.
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

        /*! Propagates sensitivities with respect to curve nodes (for
            instance, obtained by repricing a trade under bumps of the
            node values) into sensitivities with respect to the helper
            quotes of every registered curve.  The keys of the input map
            must be registered curves; each array must have one entry
            per curve node.
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
