// Portability guards for golden rules 1 and 4: identical numbers on every native
// target and on WebAssembly. Each check pins an assumption the engine relies on, so a
// target or compiler flag that breaks it fails here, not as a last-bit drift in a score.
#include <gtest/gtest.h>

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

static_assert(std::numeric_limits<double>::is_iec559, "doubles must be IEEE-754 binary64");
static_assert(CHAR_BIT == 8);
// No excess precision in intermediates (x87-style extended registers would differ).
static_assert(FLT_EVAL_METHOD == 0);

std::uint64_t bits_of(double value) {
    std::uint64_t bits = 0;
    static_assert(sizeof bits == sizeof value);
    std::memcpy(&bits, &value, sizeof bits);
    return bits;
}

// Volatile inputs stop constant folding, so the expression is evaluated by the target.
double product_minus(double a, double b, double c) {
    const volatile double va = a;
    const volatile double vb = b;
    const volatile double vc = c;
    return (va * vb) - vc;
}

TEST(Platform, MultiplyAndSubtractAreRoundedSeparately) {
    // a*a = 1 + 2^-26 + 2^-54 exactly; rounded to double it is 1 + 2^-26.
    // Rounded separately: a*a - (1 + 2^-26) == 0. Fused into one FMA (contraction): 2^-54.
    // Contraction is on by default in clang on arm64, which would make results differ
    // between Apple-silicon, x86 and wasm builds. CMake sets -ffp-contract=off.
    // Only optimized builds fuse, so CI runs this in the release preset too.
    const double a = 1.0 + 0x1p-27;
    EXPECT_EQ(product_minus(a, a, 1.0 + 0x1p-26), 0.0);
}

TEST(Platform, SqrtIsCorrectlyRounded) {
    // IEEE-754 requires correctly rounded sqrt; sqrt(2) = 0x3FF6A09E667F3BCD.
    const volatile double two = 2.0;
    EXPECT_EQ(bits_of(std::sqrt(two)), 0x3FF6A09E667F3BCDULL);
}

TEST(Platform, SizeTIsWideEnoughForEngineIndices) {
    // wasm32 has a 32-bit size_t; the engine must not assume 64 bits.
    static_assert(sizeof(std::size_t) >= 4);
    EXPECT_GE(std::numeric_limits<std::size_t>::max(), std::uint64_t{0xFFFFFFFFULL});
}

}  // namespace
