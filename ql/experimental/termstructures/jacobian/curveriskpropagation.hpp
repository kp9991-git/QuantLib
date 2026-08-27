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

/*! \file curveriskpropagation.hpp
    \brief risk propagation over curve-Jacobian blocks
*/

#ifndef quantlib_experimental_curve_risk_propagation_hpp
#define quantlib_experimental_curve_risk_propagation_hpp

#include <ql/experimental/termstructures/jacobian/curvecrossjacobian.hpp>
#include <ql/math/matrixutilities/qrdecomposition.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace QuantLib {

    namespace detail {

        //! bootstrap equations in block form
        /*! The differentiated equations are

              \f[ \sum_b A_{ab} \, dz_b = dq_a \f]

            where \f$ A_{ab} \f$ is the sensitivity of curve \f$ a \f$'s
            helper quotes to curve \f$ b \f$'s nodes. Only interacting
            off-diagonal blocks are stored.
        */
        struct CurveJacobianBlocks {
            //! cumulative free-node counts including the final total
            std::vector<Size> nodeOffset;
            //! cumulative helper counts including the final total
            std::vector<Size> quoteOffset;
            //! diagonal blocks \f$ A_{aa} \f$
            std::vector<Matrix> own;
            //! coupling blocks \f$ A_{ab} \f$ for interacting pairs
            std::map<std::pair<Size, Size>, Matrix> coupling;
            //! curves whose nodes enter each curve's helper quotes
            std::vector<std::vector<Size>> dependsOn;
            //! analytical status of each helper row
            std::vector<std::vector<bool>> analyticQuotes;

            Size size() const { return own.size(); }

            Size numNodes(Size a) const {
                return nodeOffset[a + 1] - nodeOffset[a];
            }

            Size numQuotes(Size a) const {
                return quoteOffset[a + 1] - quoteOffset[a];
            }

            //! null when the curves do not interact
            const Matrix* block(Size a, Size b) const {
                if (a == b)
                    return &own[a];
                auto i = coupling.find(std::make_pair(a, b));
                return i == coupling.end() ? nullptr : &i->second;
            }
        };

        /*! Build diagonal and coupling blocks for a curve group.
            A pair interacts when helper sensitivities reach the second curve.
            Helpers without dependency metadata are treated as reaching every
            curve and are differentiated numerically.
        */
    inline CurveJacobianBlocks curveJacobianBlocks(
        const std::vector<CurveJacobianNode>& curves,
                            const CurveCrossJacobianContext& baseContext = {}) {
            Size n = curves.size();
            CurveJacobianBlocks blocks;
            blocks.nodeOffset.assign(n + 1, 0);
            blocks.quoteOffset.assign(n + 1, 0);
            blocks.own.resize(n);
            blocks.dependsOn.resize(n);
            blocks.analyticQuotes.resize(n);

            CurveCrossJacobianContext context = baseContext;
            for (Size i = 0; i < n; ++i) {
                curves[i].ensure();
                blocks.nodeOffset[i + 1] =
                    blocks.nodeOffset[i] + curves[i].numNodes();
                blocks.quoteOffset[i + 1] =
                    blocks.quoteOffset[i] + curves[i].aliveHelpers().size();
                context.addCurve(curves[i].id, curves[i].valueDependencies);
            }

            for (Size a = 0; a < n; ++a) {
                std::set<CurveId> referenced;
                bool opaque = false;
                bool unresolved = false;
                auto helpers = curves[a].aliveHelpers();
                std::vector<QuoteSensitivities> rowSensitivities;
                rowSensitivities.reserve(helpers.size());
                for (const auto& helper : helpers) {
                    rowSensitivities.push_back(
                        helper->impliedQuoteSensitivitiesByCurve());
                    const QuoteSensitivities& s = rowSensitivities.back();
                    if (!s.available) {
                        opaque = true;
                        break;
                    }
                    for (const auto& [curve, entries] : s.sensitivities)
                        referenced.insert(curve);
                    referenced.insert(s.incomplete.begin(), s.incomplete.end());
                }
                // reuse only when every row was collected
                const std::vector<QuoteSensitivities>* rowSens =
                    rowSensitivities.size() == helpers.size()
                        ? &rowSensitivities : nullptr;
                for (const auto* target : curves[a].valueDependencies.targets())
                    referenced.insert(target);
                for (const auto* dependency : referenced)
                    unresolved = unresolved ||
                                 !context.dependencyIsResolved(dependency);

                for (Size b = 0; b < n; ++b) {
                    if (b == a)
                        continue;
                    bool interacts = opaque || unresolved;
                    for (auto i = referenced.begin();
                         !interacts && i != referenced.end(); ++i)
                        interacts = context.dependsOn(*i, curves[b].id);
                    if (interacts)
                        blocks.dependsOn[a].push_back(b);
                }

                std::vector<bool> flags;
                blocks.own[a] =
                    curveCrossJacobian(curves[a], curves[a], context, &flags,
                                       rowSens);
                blocks.analyticQuotes[a] = flags;
                for (Size b : blocks.dependsOn[a]) {
                    Matrix block =
                        curveCrossJacobian(curves[a], curves[b], context,
                                           &flags, rowSens);
                    for (Size r = 0; r < flags.size(); ++r)
                        blocks.analyticQuotes[a][r] =
                            blocks.analyticQuotes[a][r] && flags[r];
                    blocks.coupling.emplace(std::make_pair(a, b),
                                            std::move(block));
                }
            }
            return blocks;
        }

        //! node and quote risk for a curve group
        struct CurveRiskPropagation {
            //! risk from fixing one curve's nodes while the others rebootstrap
            std::vector<Array> nodeRisk;
            //! risk against helper quotes
            std::vector<Array> quoteRisk;
        };

        struct CheckedQrFactorization {
            Matrix q;
            Matrix r;
            std::vector<Size> pivot;
        };

        //! factor a square full-rank system by column-pivoted QR
        inline CheckedQrFactorization checkedQrFactorization(
                const Matrix& a,
                const char* description) {
            QL_REQUIRE(a.rows() == a.columns(),
                       description << " is not square: " << a.rows() <<
                       " rows and " << a.columns() << " columns");

            const Size n = a.rows();
            CheckedQrFactorization result{Matrix(n, n), Matrix(n, n), {}};
            if (n == 0)
                return result;

            result.pivot = qrDecomposition(a, result.q, result.r, true);
            Real largest = 0.0;
            for (Size i = 0; i < n; ++i)
                largest = std::max(largest, std::fabs(result.r[i][i]));
            QL_REQUIRE(largest > 0.0,
                       description << " is rank deficient (all QR diagonals "
                                   << "are zero)");
            const Real tolerance =
                largest * n * std::numeric_limits<Real>::epsilon();
            for (Size i = 0; i < n; ++i)
                QL_REQUIRE(std::fabs(result.r[i][i]) > tolerance,
                           description << " is rank deficient at diagonal " <<
                           i << " (|Rii|=" << std::fabs(result.r[i][i]) <<
                           ", tolerance=" << tolerance << ")");
            return result;
        }

        //! solve square full-rank systems by column-pivoted QR
        inline Matrix checkedQrSolve(const Matrix& a,
                                     const Matrix& b,
                                     const char* description) {
            QL_REQUIRE(b.rows() == a.rows(),
                       description << " right-hand side has " << b.rows() <<
                       " rows for " << a.rows() << " equations");
            CheckedQrFactorization factorization =
                checkedQrFactorization(a, description);

            const Size n = a.rows();
            if (n == 0)
                return Matrix(0, b.columns());
            Matrix y = transpose(factorization.q) * b;
            for (Size ii = n; ii > 0; --ii) {
                const Size i = ii - 1;
                for (Size k = 0; k < b.columns(); ++k) {
                    for (Size j = i + 1; j < n; ++j)
                        y[i][k] -= factorization.r[i][j] * y[j][k];
                    y[i][k] /= factorization.r[i][i];
                }
            }

            Matrix x(n, b.columns());
            for (Size i = 0; i < n; ++i)
                std::copy(y.row_begin(i), y.row_end(i),
                          x.row_begin(factorization.pivot[i]));
            return x;
        }

        //! solve one square full-rank system by column-pivoted QR
        inline Array checkedQrSolve(const Matrix& a,
                                    const Array& b,
                                    const char* description) {
            QL_REQUIRE(b.size() == a.rows(),
                       description << " right-hand side has " << b.size() <<
                       " entries for " << a.rows() << " equations");
            CheckedQrFactorization factorization =
                checkedQrFactorization(a, description);

            const Size n = a.rows();
            if (n == 0)
                return {};
            Array y = transpose(factorization.q) * b;
            for (Size ii = n; ii > 0; --ii) {
                const Size i = ii - 1;
                for (Size j = i + 1; j < n; ++j)
                    y[i] -= factorization.r[i][j] * y[j];
                y[i] /= factorization.r[i][i];
            }

            Array x(n);
            for (Size i = 0; i < n; ++i)
                x[factorization.pivot[i]] = y[i];
            return x;
        }

        //! dependency components, with dependents before their dependencies
        inline std::vector<std::vector<Size>>
        curveDependencyComponents(const CurveJacobianBlocks& blocks) {
            const Size n = blocks.size();
            std::vector<Size> order(n, 0), low(n, 0), stack;
            std::vector<bool> onStack(n, false);
            std::vector<std::vector<Size>> components;
            Size counter = 0;
            std::function<void(Size)> visit = [&](Size v) {
                order[v] = low[v] = ++counter;
                stack.push_back(v);
                onStack[v] = true;
                for (Size w : blocks.dependsOn[v]) {
                    if (order[w] == 0) {
                        visit(w);
                        low[v] = std::min(low[v], low[w]);
                    } else if (onStack[w]) {
                        low[v] = std::min(low[v], order[w]);
                    }
                }
                if (low[v] == order[v]) {
                    std::vector<Size> component;
                    Size w;
                    do {
                        w = stack.back();
                        stack.pop_back();
                        onStack[w] = false;
                        component.push_back(w);
                    } while (w != v);
                    components.push_back(std::move(component));
                }
            };
            for (Size v = 0; v < n; ++v)
                if (order[v] == 0)
                    visit(v);
            std::reverse(components.begin(), components.end());
            return components;
        }

        //! solve \f$ J^T R=S \f$ for multiple node-risk scenarios
        /*! Each matrix has nodes or quotes in rows and scenarios in columns.
            Components with zero right-hand sides are not factored.
        */
        inline std::vector<Matrix> propagateCurveNodeRiskMatrix(
                const CurveJacobianBlocks& blocks,
                const std::vector<Matrix>& directNodeRisk) {
            const Size n = blocks.size();
            QL_REQUIRE(directNodeRisk.size() == n,
                       "node risk was given for " << directNodeRisk.size() <<
                       " curves but the group holds " << n);

            Size scenarios = 0;
            if (!directNodeRisk.empty())
                scenarios = directNodeRisk.front().columns();
            std::vector<Matrix> nodeRisk = directNodeRisk;
            std::vector<Matrix> quoteRisk(n);
            for (Size i = 0; i < n; ++i) {
                QL_REQUIRE(nodeRisk[i].rows() == blocks.numNodes(i),
                           "node risk has " << nodeRisk[i].rows() <<
                           " rows for a curve with " << blocks.numNodes(i) <<
                           " nodes");
                QL_REQUIRE(nodeRisk[i].columns() == scenarios,
                           "node-risk matrices have different numbers of scenarios");
                quoteRisk[i] = Matrix(blocks.numQuotes(i), scenarios, 0.0);
            }

            std::vector<std::vector<Size>> dependents(n);
            for (Size a = 0; a < n; ++a)
                for (Size b : blocks.dependsOn[a])
                    dependents[b].push_back(a);

            std::vector<bool> solved(n, false);
            for (const auto& component : curveDependencyComponents(blocks)) {
                for (Size b : component)
                    for (Size a : dependents[b])
                        if (solved[a])
                            nodeRisk[b] -= transpose(*blocks.block(a, b)) *
                                           quoteRisk[a];

                bool active = false;
                for (Size b : component)
                    active = active || std::any_of(
                        nodeRisk[b].begin(), nodeRisk[b].end(),
                        [](Real x) { return x != 0.0; });
                if (!active) {
                    for (Size b : component)
                        solved[b] = true;
                    continue;
                }

                std::vector<Size> rowOffset(component.size() + 1, 0);
                std::vector<Size> colOffset(component.size() + 1, 0);
                for (Size i = 0; i < component.size(); ++i) {
                    rowOffset[i + 1] =
                        rowOffset[i] + blocks.numQuotes(component[i]);
                    colOffset[i + 1] =
                        colOffset[i] + blocks.numNodes(component[i]);
                }
                QL_REQUIRE(rowOffset.back() == colOffset.back(),
                           "cannot propagate curve risk through a dependency "
                           "component with " << rowOffset.back() <<
                           " helper quotes and " << colOffset.back() <<
                           " free nodes");

                Matrix jc(rowOffset.back(), colOffset.back(), 0.0);
                for (Size i = 0; i < component.size(); ++i)
                    for (Size j = 0; j < component.size(); ++j) {
                        const Matrix* A =
                            blocks.block(component[i], component[j]);
                        if (A == nullptr)
                            continue;
                        for (Size r = 0; r < A->rows(); ++r)
                            std::copy(A->row_begin(r), A->row_end(r),
                                      jc.row_begin(rowOffset[i] + r) +
                                          colOffset[j]);
                    }

                Matrix rhs(colOffset.back(), scenarios, 0.0);
                for (Size i = 0; i < component.size(); ++i)
                    for (Size r = 0; r < nodeRisk[component[i]].rows(); ++r)
                        std::copy(nodeRisk[component[i]].row_begin(r),
                                  nodeRisk[component[i]].row_end(r),
                                  rhs.row_begin(colOffset[i] + r));

                Matrix result = checkedQrSolve(
                    transpose(jc), rhs, "curve-risk Jacobian component");
                for (Size i = 0; i < component.size(); ++i)
                    for (Size r = 0; r < quoteRisk[component[i]].rows(); ++r)
                        std::copy(result.row_begin(rowOffset[i] + r),
                                  result.row_end(rowOffset[i] + r),
                                  quoteRisk[component[i]].row_begin(r));
                for (Size b : component)
                    solved[b] = true;
            }
            return quoteRisk;
        }

        //! quote-risk matrices for unit risks on one curve's nodes
        inline std::vector<Matrix> inverseCurveJacobianQuoteRisk(
                const CurveJacobianBlocks& blocks,
                Size nodeCurve) {
            QL_REQUIRE(nodeCurve < blocks.size(),
                       "curve index outside the Jacobian group");
            const Size scenarios = blocks.numNodes(nodeCurve);
            std::vector<Matrix> direct(blocks.size());
            for (Size i = 0; i < blocks.size(); ++i)
                direct[i] = Matrix(blocks.numNodes(i), scenarios, 0.0);
            for (Size i = 0; i < scenarios; ++i)
                direct[nodeCurve][i][i] = 1.0;
            return propagateCurveNodeRiskMatrix(blocks, direct);
        }

        //! selected block of the inverse full bootstrap-equation Jacobian
        inline Matrix inverseCurveJacobianBlock(
                const CurveJacobianBlocks& blocks,
                Size nodeCurve,
                Size quoteCurve) {
            QL_REQUIRE(quoteCurve < blocks.size(),
                       "curve index outside the Jacobian group");
            return transpose(
                inverseCurveJacobianQuoteRisk(blocks, nodeCurve)[quoteCurve]);
        }

        //! all quote columns for selected rows of the inverse Jacobian
        inline Matrix inverseCurveJacobianRows(
                const CurveJacobianBlocks& blocks,
                Size nodeCurve) {
            QL_REQUIRE(nodeCurve < blocks.size(),
                       "curve index outside the Jacobian group");
            const Size scenarios = blocks.numNodes(nodeCurve);
            std::vector<Matrix> quoteRisk =
                inverseCurveJacobianQuoteRisk(blocks, nodeCurve);

            Matrix result(scenarios, blocks.quoteOffset.back(), 0.0);
            for (Size b = 0; b < blocks.size(); ++b) {
                Matrix block = transpose(quoteRisk[b]);
                for (Size i = 0; i < scenarios; ++i)
                    std::copy(block.row_begin(i), block.row_end(i),
                              result.row_begin(i) + blocks.quoteOffset[b]);
            }
            return result;
        }

        /*! Solve \f$ J^T r=s \f$ by dependency component. When requested,
            zero risk fixes each curve while the others re-solve. Empty risk
            arrays are treated as zero.
        */
        inline CurveRiskPropagation propagateCurveNodeRisk(
                const CurveJacobianBlocks& blocks,
                const std::vector<Array>& directNodeRisk,
                bool computeZeroRisk) {
            Size n = blocks.size();
            QL_REQUIRE(directNodeRisk.size() == n,
                       "node risk was given for " << directNodeRisk.size() <<
                       " curves but the group holds " << n);

            CurveRiskPropagation result;
            result.nodeRisk.resize(n);
            result.quoteRisk.resize(n);
            for (Size i = 0; i < n; ++i) {
                Size nodes = blocks.numNodes(i);
                QL_REQUIRE(directNodeRisk[i].empty() ||
                               directNodeRisk[i].size() == nodes,
                           "node risk size (" << directNodeRisk[i].size() <<
                           ") does not match the number of curve nodes (" <<
                           nodes << ")");
                result.nodeRisk[i] = directNodeRisk[i].empty()
                                         ? Array(nodes, 0.0)
                                         : directNodeRisk[i];
                result.quoteRisk[i] = Array(blocks.numQuotes(i), 0.0);
            }

            std::vector<std::vector<Size>> dependents(n);
            for (Size a = 0; a < n; ++a)
                for (Size b : blocks.dependsOn[a])
                    dependents[b].push_back(a);

            std::vector<bool> solved(n, false);
            for (const auto& component : curveDependencyComponents(blocks)) {
                // risk from curves already solved
                for (Size b : component)
                    for (Size a : dependents[b])
                        if (solved[a])
                            result.nodeRisk[b] -=
                                transpose(*blocks.block(a, b)) *
                                result.quoteRisk[a];

                std::vector<Size> rowOffset(component.size() + 1, 0);
                std::vector<Size> colOffset(component.size() + 1, 0);
                for (Size i = 0; i < component.size(); ++i) {
                    rowOffset[i + 1] =
                        rowOffset[i] + blocks.numQuotes(component[i]);
                    colOffset[i + 1] =
                        colOffset[i] + blocks.numNodes(component[i]);
                }

                QL_REQUIRE(rowOffset.back() == colOffset.back(),
                           "cannot propagate curve risk through a dependency "
                           "component with " << rowOffset.back() <<
                           " helper quotes and " << colOffset.back() <<
                           " free nodes");

                Matrix jc(rowOffset.back(), colOffset.back(), 0.0);
                for (Size i = 0; i < component.size(); ++i)
                    for (Size j = 0; j < component.size(); ++j) {
                        const Matrix* A =
                            blocks.block(component[i], component[j]);
                        if (A == nullptr)
                            continue;
                        for (Size r = 0; r < A->rows(); ++r)
                            std::copy(A->row_begin(r), A->row_end(r),
                                      jc.row_begin(rowOffset[i] + r) +
                                          colOffset[j]);
                    }

                Array rhs(colOffset.back(), 0.0);
                for (Size i = 0; i < component.size(); ++i)
                    std::copy(result.nodeRisk[component[i]].begin(),
                              result.nodeRisk[component[i]].end(),
                              rhs.begin() + colOffset[i]);

                Array r = checkedQrSolve(transpose(jc), rhs,
                                         "curve-risk Jacobian component");
                for (Size i = 0; i < component.size(); ++i)
                    std::copy(r.begin() + rowOffset[i],
                              r.begin() + rowOffset[i + 1],
                              result.quoteRisk[component[i]].begin());

                // Fix each curve and re-solve the rest of the cycle.
                if (computeZeroRisk && component.size() > 1) {
                    std::vector<Array> exogenous(component.size());
                    for (Size k = 0; k < component.size(); ++k) {
                        std::vector<Size> rest;
                        for (Size i = 0; i < component.size(); ++i)
                            if (i != k)
                                rest.push_back(i);

                        std::vector<Size> rowRest(rest.size() + 1, 0);
                        std::vector<Size> colRest(rest.size() + 1, 0);
                        for (Size i = 0; i < rest.size(); ++i) {
                            rowRest[i + 1] =
                                rowRest[i] + blocks.numQuotes(component[rest[i]]);
                            colRest[i + 1] =
                                colRest[i] + blocks.numNodes(component[rest[i]]);
                        }

                        Matrix jr(rowRest.back(), colRest.back(), 0.0);
                        for (Size i = 0; i < rest.size(); ++i)
                            for (Size j = 0; j < rest.size(); ++j) {
                                const Matrix* A = blocks.block(component[rest[i]],
                                                               component[rest[j]]);
                                if (A == nullptr)
                                    continue;
                                for (Size q = 0; q < A->rows(); ++q)
                                    std::copy(A->row_begin(q), A->row_end(q),
                                              jr.row_begin(rowRest[i] + q) +
                                                  colRest[j]);
                            }

                        Array rhsRest(colRest.back(), 0.0);
                        for (Size i = 0; i < rest.size(); ++i)
                            std::copy(rhs.begin() + colOffset[rest[i]],
                                      rhs.begin() + colOffset[rest[i] + 1],
                                      rhsRest.begin() + colRest[i]);

                        Array y = checkedQrSolve(
                            transpose(jr), rhsRest,
                            "curve-risk Jacobian component without one curve");

                        Array node(blocks.numNodes(component[k]));
                        std::copy(rhs.begin() + colOffset[k],
                                  rhs.begin() + colOffset[k + 1], node.begin());
                        for (Size i = 0; i < rest.size(); ++i) {
                            const Matrix* A = blocks.block(component[rest[i]],
                                                           component[k]);
                            if (A == nullptr)
                                continue;
                            Array quote(A->rows());
                            std::copy(y.begin() + rowRest[i],
                                      y.begin() + rowRest[i + 1], quote.begin());
                            node -= transpose(*A) * quote;
                        }
                        exogenous[k] = node;
                    }
                    for (Size k = 0; k < component.size(); ++k)
                        result.nodeRisk[component[k]] = exogenous[k];
                }
                for (Size b : component)
                    solved[b] = true;
            }
            return result;
        }

    }

}

#endif
