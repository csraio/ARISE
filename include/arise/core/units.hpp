#pragma once

#include "arise/core/result.hpp"

#include <cmath>
#include <compare>
#include <string>

namespace arise::core {

/// A double tagged with a physical unit. Arithmetic is only defined within one unit,
/// and construction from a bare double is explicit, so units cannot be mixed by accident.
template <typename Tag>
class Quantity {
public:
    constexpr Quantity() noexcept = default;
    constexpr explicit Quantity(double value) noexcept : value_{value} {}

    [[nodiscard]] constexpr double value() const noexcept { return value_; }

    constexpr Quantity& operator+=(Quantity other) noexcept {
        value_ += other.value_;
        return *this;
    }
    constexpr Quantity& operator-=(Quantity other) noexcept {
        value_ -= other.value_;
        return *this;
    }

    [[nodiscard]] friend constexpr Quantity operator+(Quantity a, Quantity b) noexcept {
        return Quantity{a.value_ + b.value_};
    }
    [[nodiscard]] friend constexpr Quantity operator-(Quantity a, Quantity b) noexcept {
        return Quantity{a.value_ - b.value_};
    }
    [[nodiscard]] friend constexpr Quantity operator-(Quantity a) noexcept {
        return Quantity{-a.value_};
    }
    [[nodiscard]] friend constexpr Quantity operator*(Quantity a, double scale) noexcept {
        return Quantity{a.value_ * scale};
    }
    [[nodiscard]] friend constexpr Quantity operator*(double scale, Quantity a) noexcept {
        return Quantity{scale * a.value_};
    }
    [[nodiscard]] friend constexpr Quantity operator/(Quantity a, double divisor) noexcept {
        return Quantity{a.value_ / divisor};
    }
    /// Ratio of two quantities in the same unit is dimensionless.
    [[nodiscard]] friend constexpr double operator/(Quantity a, Quantity b) noexcept {
        return a.value_ / b.value_;
    }

    friend constexpr auto operator<=>(const Quantity&, const Quantity&) noexcept = default;

private:
    double value_{0.0};
};

struct PpmTag {};
struct HzTag {};
struct MegaHertzTag {};
struct DaltonTag {};
struct AngstromTag {};

using Ppm = Quantity<PpmTag>;              ///< Chemical shift, parts per million.
using Hz = Quantity<HzTag>;                ///< Frequency offset from the reference, hertz.
using MegaHertz = Quantity<MegaHertzTag>;  ///< Spectrometer frequency, megahertz.
using Dalton = Quantity<DaltonTag>;        ///< Mass, unified atomic mass units.
using Angstrom = Quantity<AngstromTag>;    ///< Distance, 1e-10 m.

/// True when |a - b| <= tolerance. Always false if either value is NaN.
template <typename Tag>
[[nodiscard]] constexpr bool approx_equal(Quantity<Tag> a, Quantity<Tag> b,
                                          Quantity<Tag> tolerance) noexcept {
    const double diff = a.value() - b.value();
    return (diff < 0.0 ? -diff : diff) <= tolerance.value();
}

namespace detail {

[[nodiscard]] inline Result<void> check_conversion_inputs(double value, MegaHertz spectrometer) {
    if (!std::isfinite(spectrometer.value()) || spectrometer.value() <= 0.0) {
        return Error{ErrorCode::InvalidArgument,
                     "spectrometer frequency must be finite and positive, got " +
                         std::to_string(spectrometer.value()) + " MHz"};
    }
    if (!std::isfinite(value)) {
        return Error{ErrorCode::InvalidArgument, "value to convert must be finite"};
    }
    return {};
}

}  // namespace detail

/// Frequency offset to chemical shift: delta(ppm) = delta_nu(Hz) / nu0(MHz).
[[nodiscard]] inline Result<Ppm> to_ppm(Hz offset, MegaHertz spectrometer) {
    if (auto check = detail::check_conversion_inputs(offset.value(), spectrometer); !check) {
        return std::move(check).error();
    }
    return Ppm{offset.value() / spectrometer.value()};
}

/// Chemical shift to frequency offset: delta_nu(Hz) = delta(ppm) * nu0(MHz).
[[nodiscard]] inline Result<Hz> to_hz(Ppm shift, MegaHertz spectrometer) {
    if (auto check = detail::check_conversion_inputs(shift.value(), spectrometer); !check) {
        return std::move(check).error();
    }
    return Hz{shift.value() * spectrometer.value()};
}

}  // namespace arise::core
