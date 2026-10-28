#include "arise/core/units.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <type_traits>

namespace {

using arise::core::Angstrom;
using arise::core::approx_equal;
using arise::core::Dalton;
using arise::core::ErrorCode;
using arise::core::Hz;
using arise::core::MegaHertz;
using arise::core::Ppm;
using arise::core::to_hz;
using arise::core::to_ppm;

// ---- Compile-time guarantees -------------------------------------------------

template <typename A, typename B>
concept Addable = requires(A a, B b) { a + b; };

template <typename A, typename B>
concept Comparable = requires(A a, B b) { a < b; };

// Same-unit arithmetic is allowed; mixing units or bare doubles is not.
static_assert(Addable<Ppm, Ppm>);
static_assert(!Addable<Ppm, Hz>);
static_assert(!Addable<Dalton, Angstrom>);
static_assert(!Addable<Ppm, double>);
static_assert(!Comparable<Ppm, Hz>);

// Construction from double must be explicit, and there is no implicit way back.
static_assert(std::is_constructible_v<Ppm, double>);
static_assert(!std::is_convertible_v<double, Ppm>);
static_assert(!std::is_convertible_v<Ppm, double>);
static_assert(!std::is_convertible_v<Ppm, Hz>);

// Zero-overhead wrapper: same size as double and trivially copyable.
static_assert(sizeof(Ppm) == sizeof(double));
static_assert(std::is_trivially_copyable_v<Ppm>);

// Usable in constant expressions.
static_assert((Ppm{1.5} + Ppm{2.5}).value() == 4.0);
static_assert(Ppm{} == Ppm{0.0});

// ---- Runtime behavior --------------------------------------------------------

TEST(Units, DefaultConstructsToZero) {
    EXPECT_EQ(Ppm{}.value(), 0.0);
    EXPECT_EQ(Dalton{}.value(), 0.0);
}

TEST(Units, SameUnitArithmetic) {
    const Ppm a{7.26};
    const Ppm b{1.20};
    EXPECT_DOUBLE_EQ((a + b).value(), 8.46);
    EXPECT_DOUBLE_EQ((a - b).value(), 6.06);
    EXPECT_DOUBLE_EQ((-a).value(), -7.26);
    EXPECT_DOUBLE_EQ((a * 2.0).value(), 14.52);
    EXPECT_DOUBLE_EQ((2.0 * a).value(), 14.52);
    EXPECT_DOUBLE_EQ((a / 2.0).value(), 3.63);
}

TEST(Units, CompoundAssignment) {
    Dalton mass{100.0};
    mass += Dalton{0.5};
    EXPECT_DOUBLE_EQ(mass.value(), 100.5);
    mass -= Dalton{1.0};
    EXPECT_DOUBLE_EQ(mass.value(), 99.5);
}

TEST(Units, RatioOfSameUnitIsDimensionless) {
    const double ratio = Angstrom{3.0} / Angstrom{1.5};
    EXPECT_DOUBLE_EQ(ratio, 2.0);
}

TEST(Units, ComparisonWithinUnit) {
    EXPECT_LT(Ppm{1.0}, Ppm{2.0});
    EXPECT_GT(Ppm{2.0}, Ppm{1.0});
    EXPECT_EQ(Ppm{3.0}, Ppm{3.0});
    EXPECT_NE(Ppm{3.0}, Ppm{3.1});
    EXPECT_LE(Ppm{3.0}, Ppm{3.0});
}

TEST(Units, NegativeShiftsAreValid) {
    // Shielded protons (e.g. inside a cyclophane ring) can resonate below TMS.
    const Ppm shielded{-1.5};
    EXPECT_LT(shielded, Ppm{0.0});
}

TEST(Units, ApproxEqualRespectsTolerance) {
    EXPECT_TRUE(approx_equal(Ppm{7.260}, Ppm{7.265}, Ppm{0.01}));
    EXPECT_FALSE(approx_equal(Ppm{7.260}, Ppm{7.300}, Ppm{0.01}));
    EXPECT_TRUE(approx_equal(Ppm{1.0}, Ppm{1.0}, Ppm{0.0}));
    // Symmetric in its arguments.
    EXPECT_TRUE(approx_equal(Ppm{7.265}, Ppm{7.260}, Ppm{0.01}));
}

TEST(Units, ApproxEqualIsFalseForNaN) {
    const Ppm nan{std::numeric_limits<double>::quiet_NaN()};
    EXPECT_FALSE(approx_equal(nan, nan, Ppm{1.0}));
}

// ---- Hz <-> ppm conversion ---------------------------------------------------
// delta(ppm) = delta_nu(Hz) / nu0(MHz), with delta_nu measured from the reference.

TEST(UnitsConversion, HzToPpmReferenceValue) {
    // 800 Hz offset on an 800 MHz spectrometer is exactly 1 ppm.
    const auto result = to_ppm(Hz{800.0}, MegaHertz{800.0});
    ASSERT_TRUE(result.ok());
    EXPECT_NEAR(result.value().value(), 1.0, 1e-12);
}

TEST(UnitsConversion, PpmToHzReferenceValue) {
    // CHCl3 residual at 7.26 ppm on a 400 MHz spectrometer sits 2904 Hz from TMS.
    const auto result = to_hz(Ppm{7.26}, MegaHertz{400.0});
    ASSERT_TRUE(result.ok());
    EXPECT_NEAR(result.value().value(), 2904.0, 1e-9);
}

TEST(UnitsConversion, RoundTripIsStable) {
    const Ppm original{3.1415926};
    const MegaHertz field{600.13};
    const auto hz = to_hz(original, field);
    ASSERT_TRUE(hz.ok());
    const auto back = to_ppm(hz.value(), field);
    ASSERT_TRUE(back.ok());
    EXPECT_NEAR(back.value().value(), original.value(), 1e-12);
}

TEST(UnitsConversion, NegativeOffsetsConvert) {
    const auto result = to_ppm(Hz{-400.0}, MegaHertz{400.0});
    ASSERT_TRUE(result.ok());
    EXPECT_NEAR(result.value().value(), -1.0, 1e-12);
}

TEST(UnitsConversion, ZeroFrequencyIsInvalid) {
    const auto to_p = to_ppm(Hz{100.0}, MegaHertz{0.0});
    ASSERT_FALSE(to_p.ok());
    EXPECT_EQ(to_p.error().code, ErrorCode::InvalidArgument);

    const auto to_h = to_hz(Ppm{1.0}, MegaHertz{0.0});
    ASSERT_FALSE(to_h.ok());
    EXPECT_EQ(to_h.error().code, ErrorCode::InvalidArgument);
}

TEST(UnitsConversion, NegativeFrequencyIsInvalid) {
    const auto result = to_ppm(Hz{100.0}, MegaHertz{-400.0});
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
}

TEST(UnitsConversion, NonFiniteFrequencyIsInvalid) {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(to_ppm(Hz{100.0}, MegaHertz{inf}).ok());
    EXPECT_FALSE(to_ppm(Hz{100.0}, MegaHertz{nan}).ok());
    EXPECT_FALSE(to_hz(Ppm{1.0}, MegaHertz{nan}).ok());
}

TEST(UnitsConversion, NonFiniteValueIsInvalid) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto result = to_ppm(Hz{nan}, MegaHertz{400.0});
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
    EXPECT_FALSE(to_hz(Ppm{nan}, MegaHertz{400.0}).ok());
}

TEST(UnitsConversion, ErrorMessageIsDescriptive) {
    const auto result = to_ppm(Hz{100.0}, MegaHertz{0.0});
    ASSERT_FALSE(result.ok());
    EXPECT_FALSE(result.error().message.empty());
}

}  // namespace
