#include "arise/io/correlation_csv.hpp"

#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/io/carbon_csv.hpp"
#include "arise/io/peak_csv.hpp"

#include <gtest/gtest.h>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::io::CorrelationPeak;
using arise::io::parse_hmbc_csv;
using arise::io::parse_hsqc_csv;

// Asserts failure with the given code and that the message mentions every fragment.
template <typename Parser>
void expect_error(Parser parse, std::string_view text, ErrorCode code,
                  const std::vector<std::string>& fragments) {
    const auto result = parse(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, code) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

// Illustrative shifts, not measured data. Only the shape of the data matters.

// ---- HSQC -----------------------------------------------------------------------

TEST(HsqcCsv, ParsesProtonCarbonPairsInFileOrder) {
    const auto result = parse_hsqc_csv("h_ppm,c_ppm\n7.28,128.4\n2.65,29.0\n1.24,15.6\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CorrelationPeak> expected{
        {.proton = Ppm{7.28}, .carbon = Ppm{128.4}, .intensity = std::nullopt},
        {.proton = Ppm{2.65}, .carbon = Ppm{29.0}, .intensity = std::nullopt},
        {.proton = Ppm{1.24}, .carbon = Ppm{15.6}, .intensity = std::nullopt},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(HsqcCsv, ParsesSignedIntensityFromAnEditedHsqc) {
    // Multiplicity-edited: CH2 has the opposite sign to CH and CH3.
    const auto result = parse_hsqc_csv("c_ppm,intensity,h_ppm\n128.4,1.0,7.28\n29.0,-0.8,2.65\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), 2U);
    EXPECT_EQ(result.value()[1].proton, Ppm{2.65});
    EXPECT_EQ(result.value()[1].carbon, Ppm{29.0});
    EXPECT_EQ(result.value()[1].intensity, -0.8);
}

TEST(HsqcCsv, SharesThePeakListFormat) {
    const auto result =
        parse_hsqc_csv("\xEF\xBB\xBF# HSQC, CDCl3\r\n\r\n h_ppm , c_ppm \r\n 7.28 , 128.4 \r\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), 1U);
    EXPECT_EQ(result.value()[0].carbon, Ppm{128.4});
}

TEST(HsqcCsv, BothShiftColumnsAreRequired) {
    expect_error(parse_hsqc_csv, "h_ppm\n7.28\n", ErrorCode::ParseError, {"line 1", "c_ppm"});
    expect_error(parse_hsqc_csv, "c_ppm\n128.4\n", ErrorCode::ParseError, {"line 1", "h_ppm"});
    // The 1D column name is not accepted: a 2D file must say which axis is which.
    expect_error(parse_hsqc_csv, "shift_ppm,c_ppm\n7.28,128.4\n", ErrorCode::ParseError,
                 {"line 1", "shift_ppm"});
}

TEST(HsqcCsv, AcceptsTheWindowEdgesOfBothNuclei) {
    const auto result = parse_hsqc_csv("h_ppm,c_ppm\n-10,-50\n25,300\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].proton, arise::io::kMinProtonShift);
    EXPECT_EQ(result.value()[0].carbon, arise::io::kMinCarbonShift);
    EXPECT_EQ(result.value()[1].proton, arise::io::kMaxProtonShift);
    EXPECT_EQ(result.value()[1].carbon, arise::io::kMaxCarbonShift);
}

TEST(HsqcCsv, SwappedColumnsAreCaught) {
    // 128.4 is no 1H shift: the columns were most likely exchanged.
    expect_error(parse_hsqc_csv, "h_ppm,c_ppm\n128.4,7.28\n", ErrorCode::OutOfRange,
                 {"line 2", "h_ppm", "128.4", "swapped"});
}

TEST(HsqcCsv, CarbonShiftOutsideItsWindowIsOutOfRange) {
    expect_error(parse_hsqc_csv, "h_ppm,c_ppm\n7.28,310\n", ErrorCode::OutOfRange,
                 {"line 2", "c_ppm", "310"});
}

TEST(HsqcCsv, ZeroIntensityIsOutOfRange) {
    expect_error(parse_hsqc_csv, "h_ppm,c_ppm,intensity\n7.28,128.4,0\n", ErrorCode::OutOfRange,
                 {"line 2", "intensity"});
}

TEST(HsqcCsv, HeaderWithoutPeaksIsError) {
    expect_error(parse_hsqc_csv, "h_ppm,c_ppm\n", ErrorCode::ParseError, {"no peaks"});
}

TEST(HsqcCsv, StructuralErrorsComeFromTheSharedFormat) {
    expect_error(parse_hsqc_csv, "h_ppm,c_ppm\n7.28\n", ErrorCode::ParseError,
                 {"line 2", "expected 2 fields", "got 1"});
    expect_error(parse_hsqc_csv, "h_ppm,c_ppm\n7.28,abc\n", ErrorCode::ParseError,
                 {"line 2", "c_ppm", "abc"});
}

// ---- HMBC -----------------------------------------------------------------------

TEST(HmbcCsv, ParsesLongRangeCorrelations) {
    // The CH3 protons (1.24) correlate with the CH2 (29.0) and the ipso carbon (144.2).
    const auto result = parse_hmbc_csv("h_ppm,c_ppm,intensity\n1.24,29.0,3.1\n1.24,144.2,1.2\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CorrelationPeak> expected{
        {.proton = Ppm{1.24}, .carbon = Ppm{29.0}, .intensity = 3.1},
        {.proton = Ppm{1.24}, .carbon = Ppm{144.2}, .intensity = 1.2},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(HmbcCsv, AnyFiniteIntensityIsAllowed) {
    // Magnitude-mode HMBC has no meaningful sign, and a zero is harmless.
    const auto result = parse_hmbc_csv("h_ppm,c_ppm,intensity\n1.24,29.0,0\n2.65,128.4,-0.3\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].intensity, 0.0);
    EXPECT_EQ(result.value()[1].intensity, -0.3);
}

TEST(HmbcCsv, UsesTheSameShiftWindows) {
    expect_error(parse_hmbc_csv, "h_ppm,c_ppm\n144.2,1.24\n", ErrorCode::OutOfRange,
                 {"line 2", "h_ppm", "swapped"});
    expect_error(parse_hmbc_csv, "h_ppm,c_ppm\n1.24,-55\n", ErrorCode::OutOfRange,
                 {"line 2", "c_ppm"});
}

TEST(HmbcCsv, HeaderWithoutPeaksIsError) {
    expect_error(parse_hmbc_csv, "h_ppm,c_ppm\n", ErrorCode::ParseError, {"no peaks"});
}

}  // namespace
