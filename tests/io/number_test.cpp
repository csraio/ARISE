#include "arise/io/number.hpp"

#include "arise/core/error.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string_view>

namespace {

using arise::core::ErrorCode;
using arise::io::parse_decimal;

TEST(ParseDecimal, AcceptsPlainDecimals) {
    struct Case {
        std::string_view text;
        double expected;
    };
    constexpr std::array kCases{
        Case{.text = "7.26", .expected = 7.26},     Case{.text = "0", .expected = 0.0},
        Case{.text = "-0.5", .expected = -0.5},     Case{.text = "+1", .expected = 1.0},
        Case{.text = "12", .expected = 12.0},       Case{.text = ".5", .expected = 0.5},
        Case{.text = "5.", .expected = 5.0},        Case{.text = "1e2", .expected = 100.0},
        Case{.text = "1E2", .expected = 100.0},     Case{.text = "1.5e-3", .expected = 1.5e-3},
        Case{.text = "-2.5E+1", .expected = -25.0}, Case{.text = "007", .expected = 7.0},
    };
    for (const auto& c : kCases) {
        const auto result = parse_decimal(c.text);
        ASSERT_TRUE(result.ok()) << "input: " << c.text;
        EXPECT_EQ(result.value(), c.expected) << "input: " << c.text;
    }
}

TEST(ParseDecimal, RoundsExactlyLikeTheCompiler) {
    // Exact equality on purpose: the runtime conversion must give the same bits as the
    // compiler's correctly rounded reading of the literal, on every platform.
    EXPECT_EQ(parse_decimal("0.1").value(), 0.1);
    EXPECT_EQ(parse_decimal("7.26").value(), 7.26);
    EXPECT_EQ(parse_decimal("123456.789").value(), 123456.789);
    EXPECT_EQ(parse_decimal("2.2250738585072014e-308").value(), 2.2250738585072014e-308);
    EXPECT_EQ(parse_decimal("1.7976931348623157e308").value(), 1.7976931348623157e308);
}

TEST(ParseDecimal, KeepsTheSignOfNegativeZero) {
    const auto result = parse_decimal("-0.0");
    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.value(), 0.0);
    EXPECT_TRUE(std::signbit(result.value()));
}

TEST(ParseDecimal, RejectsMalformedText) {
    constexpr std::array<std::string_view, 21> kBad{
        "",    " ",     " 1",  "1 ",  "1 2",  "+",   "-",     ".",  "e5", "1e", "1e+",
        "--1", "1.2.3", "1,5", "inf", "-inf", "nan", "0x1p3", "1f", "1d", "١",
    };
    for (const auto text : kBad) {
        const auto result = parse_decimal(text);
        ASSERT_FALSE(result.ok()) << "accepted: '" << text << "'";
        EXPECT_EQ(result.error().code, ErrorCode::ParseError) << "input: '" << text << "'";
        EXPECT_FALSE(result.error().message.empty());
    }
}

TEST(ParseDecimal, RejectsEmbeddedNul) {
    constexpr std::string_view kWithNul{
        "1\0"
        "5",
        3};
    EXPECT_FALSE(parse_decimal(kWithNul).ok());
}

TEST(ParseDecimal, RejectsOverflow) {
    const auto result = parse_decimal("1e400");
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code, ErrorCode::ParseError);
}

TEST(ParseDecimal, UnderflowRoundsToZero) {
    // Below the smallest subnormal (~4.9e-324): correctly rounded to zero, not an error.
    const auto result = parse_decimal("1e-400");
    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.value(), 0.0);
}

TEST(ParseDecimal, DoesNotReadPastTheView) {
    // The view ends before "99"; a parser that relied on NUL termination would read it.
    constexpr std::string_view kBuffer{"1.599"};
    const auto result = parse_decimal(kBuffer.substr(0, 3));
    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.value(), 1.5);
}

}  // namespace
