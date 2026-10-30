#include "arise/io/peak_csv.hpp"

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
using arise::io::kMaxProtonShift;
using arise::io::kMinProtonShift;
using arise::io::ObservedPeak;
using arise::io::parse_peak_csv;

std::vector<ObservedPeak> parse_ok(std::string_view text) {
    auto result = parse_peak_csv(text);
    EXPECT_TRUE(result.ok()) << (result.ok() ? "" : result.error().message);
    return result.ok() ? std::move(result).value() : std::vector<ObservedPeak>{};
}

// Asserts failure with the given code and that the message mentions every fragment.
void expect_error(std::string_view text, ErrorCode code,
                  const std::vector<std::string>& fragments) {
    const auto result = parse_peak_csv(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, code) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

// ---- Valid input -------------------------------------------------------------

TEST(PeakCsv, ParsesShiftsInFileOrder) {
    const auto peaks = parse_ok("shift_ppm\n7.26\n1.20\n3.85\n");
    const std::vector<ObservedPeak> expected{
        {.shift = Ppm{7.26}, .integral = std::nullopt},
        {.shift = Ppm{1.20}, .integral = std::nullopt},
        {.shift = Ppm{3.85}, .integral = std::nullopt},
    };
    EXPECT_EQ(peaks, expected);
}

TEST(PeakCsv, ParsesOptionalIntegralColumn) {
    const auto peaks = parse_ok("shift_ppm,integral\n7.26,2\n2.35,3\n");
    const std::vector<ObservedPeak> expected{
        {.shift = Ppm{7.26}, .integral = 2.0},
        {.shift = Ppm{2.35}, .integral = 3.0},
    };
    EXPECT_EQ(peaks, expected);
}

TEST(PeakCsv, ColumnsMayAppearInAnyOrder) {
    const auto peaks = parse_ok("integral,shift_ppm\n1.98,7.26\n");
    ASSERT_EQ(peaks.size(), 1U);
    EXPECT_EQ(peaks[0].shift, Ppm{7.26});
    EXPECT_EQ(peaks[0].integral, 1.98);
}

TEST(PeakCsv, ShiftsKeepExactDoubleValues) {
    // Exact equality: parsing must reproduce the compiler's value bit for bit.
    const auto peaks = parse_ok("shift_ppm\n7.26\n0.1\n");
    ASSERT_EQ(peaks.size(), 2U);
    EXPECT_EQ(peaks[0].shift.value(), 7.26);
    EXPECT_EQ(peaks[1].shift.value(), 0.1);
}

TEST(PeakCsv, AcceptsCrLfLineEndings) {
    const auto peaks = parse_ok("shift_ppm,integral\r\n7.26,1\r\n1.20,3\r\n");
    ASSERT_EQ(peaks.size(), 2U);
    EXPECT_EQ(peaks[1].shift, Ppm{1.20});
    EXPECT_EQ(peaks[1].integral, 3.0);
}

TEST(PeakCsv, AcceptsMissingFinalNewline) {
    EXPECT_EQ(parse_ok("shift_ppm\n7.26").size(), 1U);
}

TEST(PeakCsv, SkipsUtf8ByteOrderMark) {
    // Excel's "CSV UTF-8" export starts with EF BB BF.
    const auto peaks = parse_ok("\xEF\xBB\xBFshift_ppm\n7.26\n");
    ASSERT_EQ(peaks.size(), 1U);
    EXPECT_EQ(peaks[0].shift, Ppm{7.26});
}

TEST(PeakCsv, IgnoresBlankLinesAndComments) {
    const auto peaks = parse_ok(
        "# acquired in CDCl3, 400 MHz\n"
        "\n"
        "shift_ppm\n"
        "   # aromatic region\n"
        "7.26\n"
        "\n"
        "  \t\n"
        "1.20\n");
    ASSERT_EQ(peaks.size(), 2U);
    EXPECT_EQ(peaks[0].shift, Ppm{7.26});
    EXPECT_EQ(peaks[1].shift, Ppm{1.20});
}

TEST(PeakCsv, TrimsSpacesAndTabsAroundFields) {
    const auto peaks = parse_ok(" shift_ppm ,\tintegral\n 7.26\t,  2 \n");
    ASSERT_EQ(peaks.size(), 1U);
    EXPECT_EQ(peaks[0].shift, Ppm{7.26});
    EXPECT_EQ(peaks[0].integral, 2.0);
}

TEST(PeakCsv, AcceptsNegativeShifts) {
    // Internal methyls of trans-15,16-dimethyldihydropyrene, near -4.25 ppm (Boekelheide).
    const auto peaks = parse_ok("shift_ppm\n-4.25\n");
    ASSERT_EQ(peaks.size(), 1U);
    EXPECT_EQ(peaks[0].shift, Ppm{-4.25});
}

TEST(PeakCsv, AcceptsShiftsAtTheWindowEdges) {
    const auto peaks = parse_ok("shift_ppm\n-10\n25\n");
    ASSERT_EQ(peaks.size(), 2U);
    EXPECT_EQ(peaks[0].shift, kMinProtonShift);
    EXPECT_EQ(peaks[1].shift, kMaxProtonShift);
}

TEST(PeakCsv, KeepsDuplicateShifts) {
    // Two distinct protons can resonate at the same shift; nothing is merged.
    EXPECT_EQ(parse_ok("shift_ppm\n7.26\n7.26\n").size(), 2U);
}

// ---- Structural errors ----------------------------------------------------------

TEST(PeakCsv, EmptyInputIsError) {
    expect_error("", ErrorCode::ParseError, {"header"});
    expect_error("\n# only a comment\n\n", ErrorCode::ParseError, {"header"});
}

TEST(PeakCsv, HeaderWithoutPeaksIsError) {
    expect_error("shift_ppm\n", ErrorCode::ParseError, {"no peaks"});
}

TEST(PeakCsv, MissingShiftColumnIsError) {
    expect_error("integral\n2\n", ErrorCode::ParseError, {"line 1", "shift_ppm"});
}

TEST(PeakCsv, UnknownColumnIsError) {
    // A typo must not be ignored silently.
    expect_error("shift_pmm\n7.26\n", ErrorCode::ParseError, {"line 1", "shift_pmm"});
    expect_error("shift_ppm,multiplicity\n7.26,d\n", ErrorCode::ParseError,
                 {"line 1", "multiplicity"});
}

TEST(PeakCsv, ColumnNamesAreCaseSensitive) {
    expect_error("Shift_ppm\n7.26\n", ErrorCode::ParseError, {"line 1", "Shift_ppm"});
}

TEST(PeakCsv, DuplicateColumnIsError) {
    expect_error("shift_ppm,shift_ppm\n7.26,7.26\n", ErrorCode::ParseError,
                 {"line 1", "duplicate", "shift_ppm"});
}

TEST(PeakCsv, EmptyColumnNameIsError) {
    expect_error("shift_ppm,\n7.26,1\n", ErrorCode::ParseError, {"line 1"});
}

TEST(PeakCsv, RowWithTooFewFieldsIsError) {
    expect_error("shift_ppm,integral\n7.26,2\n1.20\n", ErrorCode::ParseError,
                 {"line 3", "expected 2 fields", "got 1"});
}

TEST(PeakCsv, RowWithTooManyFieldsIsError) {
    // A European decimal comma ("7,26") splits one shift into two fields.
    expect_error("shift_ppm\n7,26\n", ErrorCode::ParseError,
                 {"line 2", "expected 1 field", "got 2"});
}

TEST(PeakCsv, QuotedFieldIsError) {
    expect_error("shift_ppm\n\"7.26\"\n", ErrorCode::ParseError, {"line 2"});
}

TEST(PeakCsv, ErrorLineNumbersCountEveryPhysicalLine) {
    // Comments and blank lines still count, so the number matches the user's editor.
    expect_error("# comment\n\nshift_ppm\n7.26\n\nabc\n", ErrorCode::ParseError, {"line 6"});
    expect_error("shift_ppm\r\n7.26\r\nabc\r\n", ErrorCode::ParseError, {"line 3"});
}

// ---- Field errors ---------------------------------------------------------------

TEST(PeakCsv, NonNumericShiftIsError) {
    expect_error("shift_ppm\nabc\n", ErrorCode::ParseError, {"line 2", "shift_ppm", "abc"});
}

TEST(PeakCsv, EmptyFieldIsError) {
    expect_error("shift_ppm,integral\n7.26,\n", ErrorCode::ParseError, {"line 2", "integral"});
}

TEST(PeakCsv, NonFiniteValuesAreErrors) {
    expect_error("shift_ppm\nnan\n", ErrorCode::ParseError, {"line 2"});
    expect_error("shift_ppm\ninf\n", ErrorCode::ParseError, {"line 2"});
    expect_error("shift_ppm,integral\n7.26,1e400\n", ErrorCode::ParseError, {"line 2"});
}

// ---- Range errors ---------------------------------------------------------------

TEST(PeakCsv, CarbonShiftIsOutOfRange) {
    // 128.5 ppm is a typical aromatic 13C shift: the classic wrong-column mistake.
    expect_error("shift_ppm\n128.5\n", ErrorCode::OutOfRange, {"line 2", "128.5"});
}

TEST(PeakCsv, ShiftJustOutsideTheWindowIsOutOfRange) {
    expect_error("shift_ppm\n-10.01\n", ErrorCode::OutOfRange, {"line 2"});
    expect_error("shift_ppm\n25.01\n", ErrorCode::OutOfRange, {"line 2"});
}

TEST(PeakCsv, NonPositiveIntegralIsOutOfRange) {
    expect_error("shift_ppm,integral\n7.26,0\n", ErrorCode::OutOfRange, {"line 2", "integral"});
    expect_error("shift_ppm,integral\n7.26,-1\n", ErrorCode::OutOfRange, {"line 2", "integral"});
}

}  // namespace
