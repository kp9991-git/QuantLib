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

/*! \file curvechainrulecalculator.hpp
    \brief chain-rule propagation across curve-value dependencies
*/

#ifndef quantlib_curve_chain_rule_calculator_hpp
#define quantlib_curve_chain_rule_calculator_hpp

#include <ql/errors.hpp>
#include <ql/time/date.hpp>
#include <ql/types.hpp>
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace QuantLib {

    class TermStructure;

    namespace detail {

        using CurveId = const TermStructure*;
        using DatedCurveSensitivities = std::vector<std::pair<Date, Real>>;
        using CurveDependencyTransform =
            std::function<bool(const DatedCurveSensitivities&,
                               DatedCurveSensitivities&)>;

        class CurveChainRuleCalculator;

        //! dependencies of one curve and their sensitivity transforms
        class CurveDependencies {
          public:
            void add(CurveId target,
                     CurveDependencyTransform transform = {}) {
                QL_REQUIRE(target != nullptr, "null dependency target");
                auto existing = std::find_if(
                    entries_.begin(), entries_.end(),
                    [=](const Entry& entry) { return entry.target == target; });
                Entry entry{target, std::move(transform)};
                if (existing == entries_.end())
                    entries_.push_back(std::move(entry));
                else
                    *existing = std::move(entry);
            }

            std::vector<CurveId> targets() const {
                std::vector<CurveId> result;
                result.reserve(entries_.size());
                for (const auto& entry : entries_)
                    result.push_back(entry.target);
                return result;
            }

          private:
            struct Entry {
                CurveId target;
                //! empty when the dependency must be followed numerically
                CurveDependencyTransform transform;
            };

            std::vector<Entry> entries_;

            friend class CurveChainRuleCalculator;
        };

        //! applies the chain rule across curve-value dependencies
        class CurveChainRuleCalculator {
          public:
            void add(CurveId source,
                     CurveId target,
                     CurveDependencyTransform transform = {}) {
                QL_REQUIRE(source != nullptr, "null dependency source");
                dependencies_[source].add(target, std::move(transform));
            }

            void add(CurveId source,
                     const CurveDependencies& dependencies) {
                QL_REQUIRE(source != nullptr, "null dependency source");
                for (const auto& entry : dependencies.entries_)
                    add(source, entry.target, entry.transform);
            }

            void add(const CurveChainRuleCalculator& other) {
                for (const auto& [source, dependencies] : other.dependencies_)
                    add(source, dependencies);
            }

            std::vector<CurveId> targets(CurveId source) const {
                return outgoing(source).targets();
            }

            bool dependsOn(CurveId source, CurveId target) const {
                if (source == target)
                    return true;
                std::set<CurveId> visited;
                std::vector<CurveId> pending{source};
                while (!pending.empty()) {
                    CurveId current = pending.back();
                    pending.pop_back();
                    if (!visited.insert(current).second)
                        continue;
                    for (const auto& entry : outgoing(current).entries_) {
                        if (entry.target == target)
                            return true;
                        pending.push_back(entry.target);
                    }
                }
                return false;
            }

            bool propagate(CurveId source,
                           CurveId target,
                           const DatedCurveSensitivities& input,
                           DatedCurveSensitivities& output) const {
                std::set<CurveId> visiting;
                return propagate(source, target, input, output, visiting);
            }

          private:
            bool propagate(CurveId source,
                           CurveId target,
                           const DatedCurveSensitivities& input,
                           DatedCurveSensitivities& output,
                           std::set<CurveId>& visiting) const {
                if (source == target) {
                    output.insert(output.end(), input.begin(), input.end());
                    return true;
                }
                if (!visiting.insert(source).second)
                    return false;

                bool foundPath = false;
                DatedCurveSensitivities accumulated;
                for (const auto& entry : outgoing(source).entries_) {
                    if (entry.target != target && !dependsOn(entry.target, target))
                        continue;
                    if (!entry.transform) {
                        visiting.erase(source);
                        return false;
                    }
                    DatedCurveSensitivities direct, branch;
                    if (!entry.transform(input, direct) ||
                        !propagate(entry.target, target, direct, branch, visiting)) {
                        visiting.erase(source);
                        return false;
                    }
                    foundPath = true;
                    accumulated.insert(accumulated.end(),
                                       branch.begin(), branch.end());
                }
                visiting.erase(source);
                if (foundPath)
                    output.insert(output.end(), accumulated.begin(), accumulated.end());
                return foundPath;
            }

            const CurveDependencies& outgoing(CurveId source) const {
                static const CurveDependencies empty;
                auto i = dependencies_.find(source);
                return i == dependencies_.end() ? empty : i->second;
            }

            std::map<CurveId, CurveDependencies> dependencies_;
        };

        //! environment used while differentiating one cross-curve block
        class CurveCrossJacobianContext {
          public:
            void addCurve(
                    CurveId curve,
                    const CurveDependencies& dependencies = {}) {
                QL_REQUIRE(curve != nullptr, "null accounted curve");
                accountedCurves_.insert(curve);
                chainRule_.add(curve, dependencies);
            }

            void addDependencies(
                    CurveId curve,
                    const CurveDependencies& dependencies) {
                chainRule_.add(curve, dependencies);
            }

            void addDependencies(const CurveChainRuleCalculator& dependencies) {
                chainRule_.add(dependencies);
            }

            void addNumericallyPropagatedCurve(CurveId curve) {
                QL_REQUIRE(curve != nullptr,
                           "null numerically propagated curve");
                numericallyPropagatedCurves_.insert(curve);
            }

            void addNumericallyPropagatedCurves(
                    const std::set<CurveId>& curves) {
                for (CurveId curve : curves)
                    addNumericallyPropagatedCurve(curve);
            }

            void assumeUnlistedCurvesIndependent(bool enabled = true) {
                unlistedCurvesAreIndependent_ = enabled;
            }

            bool accountsFor(CurveId curve, CurveId differentiatedCurve) const {
                return curve == differentiatedCurve ||
                       accountedCurves_.count(curve) != 0 ||
                       (unlistedCurvesAreIndependent_ &&
                        numericallyPropagatedCurves_.count(curve) == 0);
            }

            bool dependsOn(CurveId source, CurveId target) const {
                return chainRule_.dependsOn(source, target);
            }

            bool propagate(CurveId source,
                           CurveId target,
                           const DatedCurveSensitivities& input,
                           DatedCurveSensitivities& output) const {
                return chainRule_.propagate(source, target, input, output);
            }

          private:
            std::set<CurveId> accountedCurves_;
            std::set<CurveId> numericallyPropagatedCurves_;
            bool unlistedCurvesAreIndependent_ = false;
            CurveChainRuleCalculator chainRule_;
        };

    }

}

#endif
