#include "arise/io/candidates.hpp"

#include "arise/core/error.hpp"

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::io::Candidate;
using arise::io::parse_candidate_list;

std::vector<Candidate> parse_ok(std::string_view text) {
    auto result = parse_candidate_list(text);
    EXPECT_TRUE(result.ok()) << (result.ok() ? "" : result.error().message);
    return result.ok() ? std::move(result).value() : std::vector<Candidate>{};
}

void expect_error(std::string_view text, const std::vector<std::string>& fragments) {
    const auto result = parse_candidate_list(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, ErrorCode::ParseError) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

// ---- Valid input -------------------------------------------------------------

TEST(CandidateList, ParsesSmilesAndNamesInFileOrder) {
    const auto candidates = parse_ok("CCc1ccccc1 ethylbenzene\nCC(=O)CC 2-butanone\n");
    const std::vector<Candidate> expected{
        {.smiles = "CCc1ccccc1", .name = "ethylbenzene", .line = 1},
        {.smiles = "CC(=O)CC", .name = "2-butanone", .line = 2},
    };
    EXPECT_EQ(candidates, expected);
}

TEST(CandidateList, NameIsOptional) {
    const auto candidates = parse_ok("Cc1ccccc1C\nCc1cccc(C)c1\n");
    ASSERT_EQ(candidates.size(), 2U);
    EXPECT_EQ(candidates[0].smiles, "Cc1ccccc1C");
    EXPECT_TRUE(candidates[0].name.empty());
}

TEST(CandidateList, NameIsTheRestOfTheLineTrimmed) {
    const auto candidates = parse_ok("CCO\t  ethyl alcohol (95%)  \n");
    ASSERT_EQ(candidates.size(), 1U);
    EXPECT_EQ(candidates[0].name, "ethyl alcohol (95%)");
}

TEST(CandidateList, NamesMayUseAnyUtf8) {
    const auto candidates = parse_ok("CC(C)=O Aceton (\xC3\xA4hnlich)\n");
    ASSERT_EQ(candidates.size(), 1U);
    EXPECT_EQ(candidates[0].name, "Aceton (\xC3\xA4hnlich)");
}

TEST(CandidateList, TripleBondIsDataNotAComment) {
    const auto candidates = parse_ok("# candidates from Sherlock\nC#N acetonitrile\n");
    ASSERT_EQ(candidates.size(), 1U);
    EXPECT_EQ(candidates[0].smiles, "C#N");
    EXPECT_EQ(candidates[0].line, 2U);
}

TEST(CandidateList, AcceptsTheFullSmilesAlphabet) {
    // Brackets, charges, isotopes, stereo, ring closures (incl. %nn), dots, wildcards.
    const auto candidates = parse_ok(
        "[NH4+].[Cl-]\n"
        "C[C@@H](O)CC\n"
        "F/C=C\\F\n"
        "[13CH4]\n"
        "C%10CCCCC%10\n"
        "c1ccc2ccccc2c1\n"
        "[Na+].[O-]C(=O)C\n"
        "C$C\n"
        "c:c\n"
        "*C\n");
    EXPECT_EQ(candidates.size(), 10U);
}

TEST(CandidateList, SharesTheLineRules) {
    const auto candidates = parse_ok("\xEF\xBB\xBF# list\r\n\r\n  CCO ethanol\r\n\r\nCC\r\n");
    const std::vector<Candidate> expected{
        {.smiles = "CCO", .name = "ethanol", .line = 3},
        {.smiles = "CC", .name = "", .line = 5},
    };
    EXPECT_EQ(candidates, expected);
}

TEST(CandidateList, KeepsDuplicates) {
    // Removing duplicates needs canonical SMILES, which is chem's job.
    EXPECT_EQ(parse_ok("CCO a\nCCO a\n").size(), 2U);
}

// ---- Errors -------------------------------------------------------------------

TEST(CandidateList, EmptyListIsError) {
    expect_error("", {"no candidates"});
    expect_error("# only a comment\n\n", {"no candidates"});
}

TEST(CandidateList, CommaSuggestsACsvFile) {
    expect_error("CCO ethanol\nCCO,ethanol\n", {"line 2", "','", "CSV"});
}

TEST(CandidateList, CharactersOutsideTheSmilesAlphabetAreErrors) {
    expect_error("\"CCO\" ethanol\n", {"line 1", "'\"'"});
    expect_error("CC;O\n", {"line 1", "';'"});
    expect_error("CC\xC3\xA4O\n", {"line 1", "non-ASCII"});
    // 0x80 is the first non-ASCII byte; the message must not echo the raw byte.
    expect_error(
        "C\x80"
        "C\n",
        {"line 1", "non-ASCII"});
    expect_error(
        "C\x01"
        "C\n",
        {"line 1", "control"});
}

TEST(CandidateList, HeaderLineIsError) {
    expect_error("SMILES Name\nCCO ethanol\n", {"line 1", "header"});
    expect_error("smiles\nCCO\n", {"line 1", "header"});
}

}  // namespace
