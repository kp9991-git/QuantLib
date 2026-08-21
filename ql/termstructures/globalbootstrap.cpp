/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 SoftSolutions! S.r.l.
 Copyright (C) 2025 Peter Caspers

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

#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/termstructures/globalbootstrap.hpp>

namespace QuantLib {

MultiCurveBootstrap::MultiCurveBootstrap(Real accuracy, bool analyticJacobian)
: accuracy_(accuracy), analyticJacobian_(analyticJacobian), defaultOptimizer_(true) {
    optimizer_ = ext::make_shared<LevenbergMarquardt>(accuracy, accuracy, accuracy);
    endCriteria_ = ext::make_shared<EndCriteria>(1000, 10, accuracy, accuracy, accuracy);
}

MultiCurveBootstrap::MultiCurveBootstrap(ext::shared_ptr<OptimizationMethod> optimizer,
                                         ext::shared_ptr<EndCriteria> endCriteria,
                                         bool analyticJacobian)
: optimizer_(std::move(optimizer)), endCriteria_(std::move(endCriteria)),
  analyticJacobian_(analyticJacobian) {
    constexpr auto accuracy = 1E-10;
    if (optimizer_ == nullptr) {
        optimizer_ = ext::make_shared<LevenbergMarquardt>(accuracy, accuracy, accuracy);
        // allow an analytical-Jacobian optimizer later
        defaultOptimizer_ = true;
    }
    if (endCriteria_ == nullptr)
        endCriteria_ = ext::make_shared<EndCriteria>(1000, 10, accuracy, accuracy, accuracy);
}

void MultiCurveBootstrap::add(const MultiCurveBootstrapContributor* c) {
    contributors_.push_back(c);
    c->setParentBootstrapper(shared_from_this());
}

void MultiCurveBootstrap::addObserver(Observer* o) {
    observers_.push_back(o);
}

std::set<const TermStructure*> MultiCurveBootstrap::observerTermStructures() const {
    std::set<const TermStructure*> result;
    for (auto* o : observers_)
        if (const auto* ts = dynamic_cast<const TermStructure*>(o))
            result.insert(ts);
    return result;
}

void MultiCurveBootstrap::setCostFunctionArguments(const Array& x,
                                                   const std::vector<Size>& guessSizes) const {
    std::size_t offset = 0;
    for (std::size_t c = 0; c < contributors_.size(); ++c) {
        Array tmp(guessSizes[c]);
        std::copy(std::next(x.begin(), offset), std::next(x.begin(), offset + guessSizes[c]),
                  tmp.begin());
        offset += guessSizes[c];
        contributors_[c]->setCostFunctionArgument(tmp);
    }

    for(auto *o: observers_)
        o->update();
}

// stacked cost function with optional analytical Jacobian
class MultiCurveBootstrap::StackedCostFunction : public CostFunction {
  public:
    StackedCostFunction(const MultiCurveBootstrap* b, std::vector<Size> guessSizes)
    : b_(b), guessSizes_(std::move(guessSizes)) {}
    Array values(const Array& x) const override {
        b_->setCostFunctionArguments(x, guessSizes_);

        std::vector<Array> results;
        results.reserve(b_->contributors_.size());
        for (auto& contributor : b_->contributors_) {
            results.push_back(contributor->evaluateCostFunction());
        }

        std::size_t resultSize =
            std::accumulate(results.begin(), results.end(), (std::size_t)0,
                            [](std::size_t len, const Array& a) { return len + a.size(); });

        Array result(resultSize);

        std::size_t offset = 0;
        for (auto const& r : results) {
            std::copy(r.begin(), r.end(), std::next(result.begin(), offset));
            offset += r.size();
        }

        return result;
    }
    void jacobian(Matrix& jac, const Array& x) const override {
        if (!(b_->analyticJacobian_ && b_->analyticCostJacobian(jac, x, guessSizes_)))
            CostFunction::jacobian(jac, x);
    }
  private:
    const MultiCurveBootstrap* b_;
    std::vector<Size> guessSizes_;
};

bool MultiCurveBootstrap::analyticCostJacobian(Matrix& jac,
                                               const Array& x,
                                               const std::vector<Size>& guessSizes) const {
    Size n = contributors_.size();
    setCostFunctionArguments(x, guessSizes);

    // gather contributor metadata and validate dimensions
    std::vector<detail::CurveJacobianNode> nodes(n);
    std::vector<Array> weights(n), dT(n);
    std::vector<Size> resOffset(n + 1, 0), varOffset(n + 1, 0);
    for (Size i = 0; i < n; ++i) {
        nodes[i] = contributors_[i]->jacobianNode();
        if (nodes[i].curve == nullptr)
            return false;
        weights[i] = contributors_[i]->residualWeights();
        if (weights[i].empty() || weights[i].size() != nodes[i].aliveHelpers().size())
            return false;
        Array xi(guessSizes[i]);
        std::copy(x.begin() + varOffset[i], x.begin() + varOffset[i] + guessSizes[i],
                  xi.begin());
        dT[i] = contributors_[i]->transformDerivatives(xi);
        if (dT[i].size() != guessSizes[i] || nodes[i].numNodes() != guessSizes[i])
            return false;
        resOffset[i + 1] = resOffset[i] + weights[i].size();
        varOffset[i + 1] = varOffset[i] + guessSizes[i];
    }
    if (jac.rows() != resOffset[n] || jac.columns() != varOffset[n])
        return false;

    // known dependent wrappers require numerical propagation
    detail::CurveCrossJacobianContext context;
    context.addNumericallyPropagatedCurves(observerTermStructures());
    context.assumeUnlistedCurvesIndependent();
    for (const auto& node : nodes)
        context.addCurve(node.id, node.valueDependencies);

    // differentiate w_i * (quote_i - impliedQuote_i) through curve nodes
    for (Size i = 0; i < n; ++i) {
        for (Size j = 0; j < n; ++j) {
            std::vector<bool> analytic;
            Matrix Jq = detail::curveCrossJacobian(nodes[i], nodes[j],
                                                   context, &analytic);
            for (auto flag : analytic)  // NOLINT(readability-use-anyofallof)
                if (!flag)
                    return false;
            for (Size r = 0; r < Jq.rows(); ++r)
                for (Size c = 0; c < Jq.columns(); ++c)
                    jac[resOffset[i] + r][varOffset[j] + c] =
                        -weights[i][r] * Jq[r][c] * dT[j][c];
        }
    }
    return true;
}

void MultiCurveBootstrap::runMultiCurveBootstrap() {

    std::vector<Size> guessSizes;
    std::vector<Real> globalGuess;

    for (auto const& c : contributors_) {
        Array guess = c->setupCostFunction();
        globalGuess.insert(globalGuess.end(), guess.begin(), guess.end());
        guessSizes.push_back(guess.size());
    }

    StackedCostFunction costFunction(this, guessSizes);
    Array guess(globalGuess.begin(), globalGuess.end());

    ext::shared_ptr<OptimizationMethod> optimizer = optimizer_;
    if (analyticJacobian_ && !guess.empty()) {
        Matrix probe(costFunction.values(guess).size(), guess.size(), 0.0);
        QL_REQUIRE(analyticCostJacobian(probe, guess, guessSizes),
                   "the analytical Jacobian was requested, but it is not "
                   "available for this multi-curve bootstrap");
        if (defaultOptimizer_)
            optimizer = ext::make_shared<LevenbergMarquardt>(
                accuracy_, accuracy_, accuracy_, /*useCostFunctionsJacobian=*/true);
        QL_REQUIRE(optimizer->usesCostFunctionJacobian(),
                   "the analytical Jacobian was requested, but the supplied "
                   "optimizer does not consume CostFunction::jacobian()");
    }

    NoConstraint noConstraint;
    Problem problem(costFunction, noConstraint, guess);
    EndCriteria::Type endType = optimizer->minimize(problem, *endCriteria_);

    QL_REQUIRE(
        EndCriteria::succeeded(endType),
        "global bootstrap failed to minimize to required accuracy (during multi curve bootstrap): "
            << endType);

    // set all contributors to valid

    for (auto const& c : contributors_)
        c->setToValid();
}

} // namespace QuantLib
