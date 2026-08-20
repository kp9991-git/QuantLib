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
#include <algorithm>
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
                    bool firstTied = false;
                    if constexpr (hasFirstDataPointFlag<Traits>)
                        firstTied = Traits::firstDataPointTracksSecond;

                    Time tMax = times.back();
                    std::vector<Real> row(cols);
                    for (Size i = 0; i < rows; ++i) {
                        auto sensitivities = alive[i]->impliedQuoteSensitivities();
                        if (sensitivities.empty())
                            continue;
                        std::fill(row.begin(), row.end(), 0.0);
                        bool ok = true;
                        for (const auto& [t, dQdP] : sensitivities) {
                            // beyond the last node the curve might extrapolate
                            // differently than the interpolation does
                            if (t > tMax) {
                                ok = false;
                                break;
                            }
                            auto weights = interpolation.nodeWeights(t, true);
                            if (weights.empty()) {
                                ok = false;
                                break;
                            }
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
                        if (ok) {
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

                    // restore element-wise: assigning the vector would
                    // invalidate the iterators stored in the interpolation
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

    }

}

#endif
