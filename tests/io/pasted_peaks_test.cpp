#include "arise/io/pasted_peaks.hpp"

#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/io/carbon_csv.hpp"

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
using arise::io::parse_pasted_carbon;
using arise::io::parse_pasted_dept;

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

// Illustrative shifts, not measured data: only the text layout matters here.

// ---- 13C ------------------------------------------------------------------------

TEST(PastedCarbon, OneShiftPerLineInPastedOrder) {
    const auto result = parse_pasted_carbon("170.5\n29\n128.25\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CarbonPeak> expected{
        {.shift = Ppm{170.5}, .intensity = std::nullopt},
        {.shift = Ppm{29.0}, .intensity = std::nullopt},
        {.shift = Ppm{128.25}, .intensity = std::nullopt},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(PastedCarbon, ReadsAnIntensityAfterTabsSemicolonsOrSpaces) {
    const auto result = parse_pasted_carbon("170.5\t1200\n29;-3.5\n128.25   40\n77 \t 1e3\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CarbonPeak> expected{
        {.shift = Ppm{170.5}, .intensity = 1200.0},
        {.shift = Ppm{29.0}, .intensity = -3.5},
        {.shift = Ppm{128.25}, .intensity = 40.0},
        {.shift = Ppm{77.0}, .intensity = 1000.0},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(PastedCarbon, IntensityIsOptionalLineByLine) {
    const auto result = parse_pasted_carbon("170.5\t1200\n29\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<CarbonPeak> expected{
        {.shift = Ppm{170.5}, .intensity = 1200.0},
        {.shift = Ppm{29.0}, .intensity = std::nullopt},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(PastedCarbon, GivesTheSamePeaksAsTheCsvReader) {
    const auto pasted = parse_pasted_carbon("170.5\t1200\n29\t-3.5\n128.25\t40\n");
    const auto csv = parse_carbon_csv("shift_ppm,intensity\n170.5,1200\n29,-3.5\n128.25,40\n");
    ASSERT_TRUE(pasted.ok()) << pasted.error().message;
    ASSERT_TRUE(csv.ok()) << csv.error().message;
    EXPECT_EQ(pasted.value(), csv.value());
}

TEST(PastedCarbon, SkipsBlankLinesCommentsAndSurroundingSpace) {
    const auto result = parse_pasted_carbon("\n  # peak table\n\n  170.5  \r\n\t29\t\r\n\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    ASSERT_EQ(result.value().size(), 2U);
    EXPECT_EQ(result.value()[0].shift, Ppm{170.5});
    EXPECT_EQ(result.value()[1].shift, Ppm{29.0});
}

TEST(PastedCarbon, LineNumbersCountEachLineEndOnce) {
    // CRLF is one line end, as are a lone LF and a lone CR.
    expect_error(parse_pasted_carbon, "29\r\n30\r\nabc\r\n", ErrorCode::ParseError, {"line 3"});
    expect_error(parse_pasted_carbon, "29\r30\rabc\r", ErrorCode::ParseError, {"line 3"});
    expect_error(parse_pasted_carbon, "29\n\r\nabc\n", ErrorCode::ParseError, {"line 3"});
}

TEST(PastedCarbon, AcceptsOldMacLineEnds) {
    const auto result = parse_pasted_carbon("170.5\r29\r");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value().size(), 2U);
}

TEST(PastedCarbon, SkipsAHeaderOnTheFirstLineOnly) {
    const auto result = parse_pasted_carbon("ppm\tIntensity\n170.5\t1200\n29\t40\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value().size(), 2U);
    // A header after blank lines and comments is still the first line that counts.
    EXPECT_TRUE(parse_pasted_carbon("\n# exported\nShift (ppm)\n29\n").ok());
    // A second text line is not a header.
    expect_error(parse_pasted_carbon, "ppm\n170.5\nppm\n", ErrorCode::ParseError, {"line 3"});
}

TEST(PastedCarbon, AFirstLineWithAnyNumberIsDataNotAHeader) {
    // "170.5 abc" has a number, so it is read as data and its second field fails.
    expect_error(parse_pasted_carbon, "170.5\tabc\n29\n", ErrorCode::ParseError, {"line 1"});
}

TEST(PastedCarbon, CommaDecimalsAreAnErrorNotTwoNumbers) {
    expect_error(parse_pasted_carbon, "170.5\n29,5\n", ErrorCode::ParseError, {"line 2"});
}

TEST(PastedCarbon, AMistypedFirstPeakIsAnErrorNotAHeader) {
    // Skipping it as a header would silently lose the peak.
    expect_error(parse_pasted_carbon, "29,5\n170.5\n", ErrorCode::ParseError, {"line 1"});
    expect_error(parse_pasted_carbon, "-29x\n170.5\n", ErrorCode::ParseError, {"line 1"});
    expect_error(parse_pasted_carbon, ".5ppm\n170.5\n", ErrorCode::ParseError, {"line 1"});
    // A first line without a letter is not a header either.
    expect_error(parse_pasted_carbon, "***\n170.5\n", ErrorCode::ParseError, {"line 1"});
    expect_error(parse_pasted_carbon, "---\n170.5\n", ErrorCode::ParseError, {"line 1"});
}

TEST(PastedCarbon, AThirdColumnIsAnError) {
    // e.g. a peak-number column pasted along: refusing it beats guessing the columns.
    expect_error(parse_pasted_carbon, "1\t170.5\t1200\n", ErrorCode::ParseError, {"line 1"});
}

TEST(PastedCarbon, NoPeaksIsAnError) {
    expect_error(parse_pasted_carbon, "", ErrorCode::ParseError, {});
    expect_error(parse_pasted_carbon, "\n\n# nothing\n", ErrorCode::ParseError, {});
    expect_error(parse_pasted_carbon, "ppm\tIntensity\n", ErrorCode::ParseError, {});
}

TEST(PastedCarbon, ShiftLimitsAreTheCsvReadersLimits) {
    const auto edges = parse_pasted_carbon("-50\n300\n");
    ASSERT_TRUE(edges.ok()) << edges.error().message;
    EXPECT_EQ(edges.value()[0].shift, kMinCarbonShift);
    EXPECT_EQ(edges.value()[1].shift, kMaxCarbonShift);
    expect_error(parse_pasted_carbon, "29\n300.5\n", ErrorCode::OutOfRange, {"line 2"});
    expect_error(parse_pasted_carbon, "-50.5\n", ErrorCode::OutOfRange, {"line 1"});
}

TEST(PastedCarbon, NonNumbersAreParseErrors) {
    expect_error(parse_pasted_carbon, "29\nnan\n", ErrorCode::ParseError, {"line 2"});
    expect_error(parse_pasted_carbon, "29\n0x1A\n", ErrorCode::ParseError, {"line 2"});
    expect_error(parse_pasted_carbon, "29\n29 ppm\n", ErrorCode::ParseError, {"line 2"});
}

// ---- DEPT -----------------------------------------------------------------------

TEST(PastedDept, ReadsSignedIntensities) {
    const auto result = parse_pasted_dept("29\t-850\n128.25 1200\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<DeptPeak> expected{
        {.shift = Ppm{29.0}, .intensity = -850.0},
        {.shift = Ppm{128.25}, .intensity = 1200.0},
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(PastedDept, GivesTheSamePeaksAsTheCsvReader) {
    const auto pasted = parse_pasted_dept("29\t-850\n128.25\t1200\n");
    const auto csv = parse_dept_csv("shift_ppm,intensity\n29,-850\n128.25,1200\n");
    ASSERT_TRUE(pasted.ok()) << pasted.error().message;
    ASSERT_TRUE(csv.ok()) << csv.error().message;
    EXPECT_EQ(pasted.value(), csv.value());
}

TEST(PastedDept, EveryLineNeedsAnIntensity) {
    expect_error(parse_pasted_dept, "29\t-850\n128.25\n", ErrorCode::ParseError, {"line 2"});
}

TEST(PastedDept, ZeroIntensityIsOutOfRange) {
    // A zero has no phase, as in parse_dept_csv.
    expect_error(parse_pasted_dept, "29\t0\n", ErrorCode::OutOfRange, {"line 1"});
    expect_error(parse_pasted_dept, "29\t-0.0\n", ErrorCode::OutOfRange, {"line 1"});
}

TEST(PastedDept, ShiftLimitsAreTheCsvReadersLimits) {
    expect_error(parse_pasted_dept, "310\t1\n", ErrorCode::OutOfRange, {"line 1"});
}

TEST(PastedDept, NoPeaksIsAnEmptySpectrum) {
    for (const std::string_view text : {"", "\n# none\n", "ppm\tIntensity\n"}) {
        const auto result = parse_pasted_dept(text);
        ASSERT_TRUE(result.ok()) << result.error().message;
        EXPECT_TRUE(result.value().empty());
    }
}

TEST(PastedDept, AThirdColumnIsAnError) {
    expect_error(parse_pasted_dept, "29\t-850\t3\n", ErrorCode::ParseError, {"line 1"});
}

}  // namespace
