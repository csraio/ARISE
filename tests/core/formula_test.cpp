#include "arise/core/formula.hpp"

#include "arise/core/error.hpp"

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Formula;

Formula parse_ok(std::string_view text) {
    auto result = Formula::parse(text);
    EXPECT_TRUE(result.ok()) << text << ": " << (result.ok() ? "" : result.error().message);
    return result.ok() ? std::move(result).value() : Formula{};
}

void expect_error(std::string_view text, ErrorCode code,
                  const std::vector<std::string>& fragments) {
    const auto result = Formula::parse(text);
    ASSERT_FALSE(result.ok()) << "accepted '" << text << "'";
    EXPECT_EQ(result.error().code, code) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

// ---- Parsing -------------------------------------------------------------------

TEST(Formula, ParsesHillNotation) {
    const auto formula = parse_ok("C4H8O2");
    EXPECT_EQ(formula.count("C"), 4);
    EXPECT_EQ(formula.count("H"), 8);
    EXPECT_EQ(formula.count("O"), 2);
    EXPECT_EQ(formula.count("N"), 0);
    EXPECT_EQ(formula.hill(), "C4H8O2");
}

TEST(Formula, MissingCountMeansOne) {
    const auto formula = parse_ok("CH4O");
    EXPECT_EQ(formula.count("C"), 1);
    EXPECT_EQ(formula.count("O"), 1);
    EXPECT_EQ(formula.hill(), "CH4O");
}

TEST(Formula, AnyOrderNormalisesToHill) {
    EXPECT_EQ(parse_ok("O2H8C4").hill(), "C4H8O2");
    EXPECT_EQ(parse_ok("NC5H5").hill(), "C5H5N");
    EXPECT_EQ(parse_ok("SOC2H6").hill(), "C2H6OS");
}

TEST(Formula, WithoutCarbonEverythingIsAlphabetical) {
    // Hill convention: hydrogen chloride is ClH, water is H2O.
    EXPECT_EQ(parse_ok("HCl").hill(), "ClH");
    EXPECT_EQ(parse_ok("OH2").hill(), "H2O");
}

TEST(Formula, CondensedFormulasSumRepeatedElements) {
    EXPECT_EQ(parse_ok("CH3COOH").hill(), "C2H4O2");
    EXPECT_EQ(parse_ok("CH3CH2OCH2CH3").hill(), "C4H10O");
}

TEST(Formula, ReadsTwoLetterSymbolsAndMultiDigitCounts) {
    const auto formula = parse_ok("CH2Cl2");
    EXPECT_EQ(formula.count("Cl"), 2);
    EXPECT_EQ(formula.count("C"), 1);
    EXPECT_EQ(formula.hill(), "CH2Cl2");
    EXPECT_EQ(parse_ok("C15H24O").count("C"), 15);
    EXPECT_EQ(parse_ok("C6H5Br").hill(), "C6H5Br");
    EXPECT_EQ(parse_ok("C8H20Si").count("Si"), 1);
}

TEST(Formula, IgnoresSurroundingWhitespace) {
    EXPECT_EQ(parse_ok("  C4H8O2\n").hill(), "C4H8O2");
}

TEST(Formula, EqualFormulasCompareEqual) {
    EXPECT_EQ(parse_ok("C4H8O2"), parse_ok("O2C4H8"));
    EXPECT_NE(parse_ok("C4H8O2"), parse_ok("C4H8O"));
}

// ---- Parse errors --------------------------------------------------------------

TEST(Formula, EmptyInputIsError) {
    expect_error("", ErrorCode::ParseError, {"empty"});
    expect_error("   ", ErrorCode::ParseError, {"empty"});
}

TEST(Formula, LowercaseStartIsError) {
    expect_error("c4h8o2", ErrorCode::ParseError, {"capital"});
}

TEST(Formula, UnknownElementIsError) {
    expect_error("C4Xx2", ErrorCode::ParseError, {"Xx"});
    expect_error("C2H6Hg", ErrorCode::ParseError, {"Hg"});
}

TEST(Formula, CobaltProbablyMeantCarbonAndOxygen) {
    // "Co" is cobalt; a chemist typing CO2 in a hurry gets a hint, not a silent misread.
    expect_error("Co2", ErrorCode::ParseError, {"Co", "CO"});
}

TEST(Formula, BadCountsAreErrors) {
    expect_error("C0H4", ErrorCode::ParseError, {"0"});
    expect_error("C04H8", ErrorCode::ParseError, {"04"});
    expect_error("2C", ErrorCode::ParseError, {"2"});
    expect_error("C1000", ErrorCode::OutOfRange, {"1000"});
}

TEST(Formula, UnsupportedNotationIsError) {
    expect_error("C(CH3)4", ErrorCode::ParseError, {"'('"});
    expect_error("C4H8O2+", ErrorCode::ParseError, {"'+'"});
    expect_error("C4 H8 O2", ErrorCode::ParseError, {"' '"});
}

// ---- Degrees of unsaturation ---------------------------------------------------

TEST(Formula, DegreesOfUnsaturationFromStandardValences) {
    EXPECT_EQ(parse_ok("C4H8O2").degrees_of_unsaturation().value(), 1);  // one C=O
    EXPECT_EQ(parse_ok("C6H6").degrees_of_unsaturation().value(), 4);    // ring + 3 C=C
    EXPECT_EQ(parse_ok("C5H5N").degrees_of_unsaturation().value(), 4);   // pyridine
    EXPECT_EQ(parse_ok("C2H6O").degrees_of_unsaturation().value(), 0);   // ethanol
    EXPECT_EQ(parse_ok("CH2Cl2").degrees_of_unsaturation().value(), 0);
    EXPECT_EQ(parse_ok("C8H8O2").degrees_of_unsaturation().value(), 5);  // e.g. methyl benzoate
    EXPECT_EQ(parse_ok("C15H24O").degrees_of_unsaturation().value(), 4);
}

TEST(Formula, OddElectronCountIsInvalid) {
    const auto methyl = parse_ok("CH3").degrees_of_unsaturation();
    ASSERT_FALSE(methyl.ok());
    EXPECT_EQ(methyl.error().code, ErrorCode::InvalidArgument);
    EXPECT_NE(methyl.error().message.find("CH3"), std::string::npos);
}

TEST(Formula, TooManyHydrogensIsInvalid) {
    const auto result = parse_ok("CH6").degrees_of_unsaturation();
    ASSERT_FALSE(result.ok());
    EXPECT_EQ(result.error().code, ErrorCode::InvalidArgument);
}

}  // namespace
