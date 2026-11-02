#include "arise/io/carbon_csv.hpp"

#include "arise/core/error.hpp"
#include "arise/core/units.hpp"

#include <gtest/gtest.h>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::io::CarbonPeak;
using arise::io::DeptPeak;
using arise::io::kMaxCarbonShift;
using arise::io::kMinCarbonShift;
using arise::io::parse_carbon_csv;
using arise::io::parse_dept_csv;

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

// ---- 13C ------------------------------------------------------------------------
// Ethylbenzene in CDCl3, approximately (SDBS): ipso C 144.2, aromatic CH 128.4, 127.9,
// 125.7, CH2 29.0, CH3 15.6 ppm. Only the shape of the data matters here.

TEST(CarbonCsv, ParsesShiftsInFileOrder) {
    const auto result = parse_carbon_csv("shift_ppm\n144.2\n128.4\n29.0\n15.6\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CarbonPeak> expected{
        {.shift = Ppm{144.2}, .intensity = std::nullopt},
        {.shift = Ppm{128.4}, .intensity = std::nullopt},
        {.shift = Ppm{29.0}, .intensity = std::nullopt},
        {.shift = Ppm{15.6}, .intensity = std::nullopt},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(CarbonCsv, ParsesOptionalIntensityColumn) {
    const auto result = parse_carbon_csv("intensity,shift_ppm\n0.31,144.2\n1.0,128.4\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), 2U);
    EXPECT_EQ(result.value()[0].shift, Ppm{144.2});
    EXPECT_EQ(result.value()[0].intensity, 0.31);
}

TEST(CarbonCsv, SharesThePeakListFormat) {
    // Same conventions as every peak list: BOM, CRLF, comments, blank lines, trimming.
    const auto result =
        parse_carbon_csv("\xEF\xBB\xBF# 13C{1H}, CDCl3\r\n\r\n shift_ppm \r\n 77.16 \r\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), 1U);
    EXPECT_EQ(result.value()[0].shift, Ppm{77.16});
}

TEST(CarbonCsv, AcceptsShiftsAtTheWindowEdges) {
    const auto result = parse_carbon_csv("shift_ppm\n-50\n300\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].shift, kMinCarbonShift);
    EXPECT_EQ(result.value()[1].shift, kMaxCarbonShift);
}

TEST(CarbonCsv, ShiftOutsideTheWindowIsOutOfRange) {
    expect_error(parse_carbon_csv, "shift_ppm\n300.5\n", ErrorCode::OutOfRange,
                 {"line 2", "300.5", "13C"});
    expect_error(parse_carbon_csv, "shift_ppm\n-50.5\n", ErrorCode::OutOfRange, {"line 2"});
    // An offset in Hz at 100 MHz, the classic wrong-unit mistake.
    expect_error(parse_carbon_csv, "shift_ppm\n12840\n", ErrorCode::OutOfRange, {"line 2"});
}

TEST(CarbonCsv, StructuralErrorsComeFromTheSharedFormat) {
    expect_error(parse_carbon_csv, "", ErrorCode::ParseError, {"header"});
    expect_error(parse_carbon_csv, "shift_ppm\n", ErrorCode::ParseError, {"no peaks"});
    expect_error(parse_carbon_csv, "shift_ppm,integral\n128.4,1\n", ErrorCode::ParseError,
                 {"line 1", "integral"});
    expect_error(parse_carbon_csv, "shift_ppm\nabc\n", ErrorCode::ParseError, {"line 2", "abc"});
}

TEST(CarbonCsv, ZeroAndNegativeIntensitiesAreAllowed) {
    // 13C intensities carry no meaning ARISE relies on; any finite value is kept.
    const auto result = parse_carbon_csv("shift_ppm,intensity\n128.4,0\n29.0,-0.2\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].intensity, 0.0);
    EXPECT_EQ(result.value()[1].intensity, -0.2);
}

// ---- DEPT -----------------------------------------------------------------------

TEST(DeptCsv, ParsesSignedIntensities) {
    // DEPT-135 of ethylbenzene: aromatic CH and CH3 up, CH2 down, ipso C absent.
    const auto result =
        parse_dept_csv("shift_ppm,intensity\n128.4,1.0\n127.9,0.9\n29.0,-0.8\n15.6,0.7\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<DeptPeak> expected{
        {.shift = Ppm{128.4}, .intensity = 1.0},
        {.shift = Ppm{127.9}, .intensity = 0.9},
        {.shift = Ppm{29.0}, .intensity = -0.8},
        {.shift = Ppm{15.6}, .intensity = 0.7},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(DeptCsv, IntensityColumnIsRequired) {
    // Without the sign a DEPT-135 peak list cannot tell CH2 from CH or CH3.
    expect_error(parse_dept_csv, "shift_ppm\n29.0\n", ErrorCode::ParseError,
                 {"line 1", "intensity"});
}

TEST(DeptCsv, ZeroIntensityIsOutOfRange) {
    // A zero has no phase, so the peak cannot be classified.
    expect_error(parse_dept_csv, "shift_ppm,intensity\n29.0,0\n", ErrorCode::OutOfRange,
                 {"line 2", "intensity"});
    expect_error(parse_dept_csv, "shift_ppm,intensity\n29.0,-0.0\n", ErrorCode::OutOfRange,
                 {"line 2", "intensity"});
}

TEST(DeptCsv, ShiftOutsideTheCarbonWindowIsOutOfRange) {
    expect_error(parse_dept_csv, "shift_ppm,intensity\n310,1\n", ErrorCode::OutOfRange,
                 {"line 2", "310"});
}

TEST(DeptCsv, HeaderOnlyIsAnEmptySpectrum) {
    // A DEPT-90 of a molecule without CH carbons (e.g. 2-butanone) has no peaks at all.
    const auto result = parse_dept_csv("# DEPT-90\nshift_ppm,intensity\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_TRUE(result.value().empty());
}

TEST(DeptCsv, StructuralErrorsComeFromTheSharedFormat) {
    expect_error(parse_dept_csv, "", ErrorCode::ParseError, {"header"});
    expect_error(parse_dept_csv, "shift_ppm,intensity\n29.0\n", ErrorCode::ParseError,
                 {"line 2", "expected 2 fields", "got 1"});
    expect_error(parse_dept_csv, "shift_ppm,intensity\n29.0,nan\n", ErrorCode::ParseError,
                 {"line 2"});
}

}  // namespace
