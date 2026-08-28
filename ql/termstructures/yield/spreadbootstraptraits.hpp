/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#ifndef quantlib_spread_bootstrap_traits_hpp
#define quantlib_spread_bootstrap_traits_hpp

#include <ql/termstructures/yield/bootstraptraits.hpp>
#include <ql/termstructures/yield/spreaddiscountcurve.hpp>

namespace QuantLib::detail {

    template <class Traits>
    struct SpreadTraits;

    //! Spread Discount-curve traits
    template <>
    struct SpreadTraits<Discount> : Discount {
        // interpolated curve type
        template <class Interpolator>
        struct curve {
            typedef InterpolatedSpreadDiscountCurve<Interpolator> type;
        };

        template <class C>
        static Real discountFactorDerivative(Time t, const C* c) {
            // The curve data are multiplicative discount spreads, while
            // helper sensitivities are with respect to the final discount.
            return c->baseCurve()->discount(t);
        }

        template <class C>
        static std::vector<std::pair<Size, Real>> extrapolationNodeWeights(
            Time t, const C* c,
                                 const Interpolation& interpolation) {
            // The inherited weights scale with the full discount
            // B(t)*S(t), but the curve nodes only drive the spread
            // factor S(t). discountFactorDerivative() supplies B(t),
            // so divide the base discount back out.
            auto weights =
                Discount::extrapolationNodeWeights(t, c, interpolation);
            DiscountFactor baseDiscount = c->baseCurve()->discount(t, true);
            for (auto& [j, w] : weights)
                w /= baseDiscount;
            return weights;
        }

        template <class C, class Input, class Output>
        static bool transformBaseCurveSensitivities(
                const C* c, const Input& input, Output& output) {
            output.reserve(output.size() + input.size());
            for (const auto& [date, sensitivity] : input) {
                DiscountFactor baseDiscount =
                    c->baseCurve()->discount(date, true);
                DiscountFactor spread = c->discount(date, true) / baseDiscount;
                output.emplace_back(date, sensitivity * spread);
            }
            return true;
        }
    };

}

#endif
