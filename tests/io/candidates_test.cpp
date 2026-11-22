#include "arise/io/candidates.hpp"

#include "arise/core/error.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::io::Candidate;
using arise::io::parse_candidate_list;
using arise::test::compound_smiles;
using arise::test::syntax_case;

// `smiles` with `inserted` placed before its character at `position`.
std::string with_inserted(std::string smiles, std::size_t position, std::string_view inserted) {
    return smiles.insert(position, inserted);
}

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
    const auto ethylbenzene = syntax_case("aromatic lowercase ethylbenzene");
    const auto butanone = compound_smiles("Methyl Ethyl Ketone");
    const auto candidates = parse_ok(ethylbenzene + " ethylbenzene\n" + butanone + " 2-butanone\n");
    const std::vector<Candidate> expected{
        {.smiles = ethylbenzene, .name = "ethylbenzene", .line = 1},
        {.smiles = butanone, .name = "2-butanone", .line = 2},
    };
    EXPECT_EQ(candidates, expected);
}

TEST(CandidateList, NameIsOptional) {
    const auto o_xylene = syntax_case("aromatic lowercase o-xylene");
    const auto candidates =
        parse_ok(o_xylene + "\n" + syntax_case("aromatic lowercase m-xylene") + "\n");
    ASSERT_EQ(candidates.size(), 2U);
    EXPECT_EQ(candidates[0].smiles, o_xylene);
    EXPECT_TRUE(candidates[0].name.empty());
}

TEST(CandidateList, NameIsTheRestOfTheLineTrimmed) {
    const auto candidates = parse_ok(compound_smiles("Ethanol") + "\t  ethyl alcohol (95%)  \n");
    ASSERT_EQ(candidates.size(), 1U);
    EXPECT_EQ(candidates[0].name, "ethyl alcohol (95%)");
}

TEST(CandidateList, NamesMayUseAnyUtf8) {
    const auto candidates = parse_ok(compound_smiles("Acetone") + " Aceton (\xC3\xA4hnlich)\n");
    ASSERT_EQ(candidates.size(), 1U);
    EXPECT_EQ(candidates[0].name, "Aceton (\xC3\xA4hnlich)");
}

TEST(CandidateList, TripleBondIsDataNotAComment) {
    // Hydrogen cyanide's SMILES contains '#', which must not start a comment.
    const auto cyanide = compound_smiles("Hydrogen Cyanide");
    ASSERT_NE(cyanide.find('#'), std::string::npos);
    const auto candidates =
        parse_ok("# candidates from Sherlock\n" + cyanide + " hydrogen cyanide\n");
    ASSERT_EQ(candidates.size(), 1U);
    EXPECT_EQ(candidates[0].smiles, cyanide);
    EXPECT_EQ(candidates[0].line, 2U);
}

TEST(CandidateList, AcceptsTheFullSmilesAlphabet) {
    // Brackets, charges, isotopes, stereo, ring closures (incl. %nn), dots, wildcards.
    std::string text;
    for (const auto& smiles : {
             compound_smiles("Ammonium Chloride"),
             compound_smiles("(-)-2-Butanol"),
             compound_smiles("cis-1,2-Difluoroethylene"),
             compound_smiles("Carbon-13C"),
             syntax_case("two-digit ring closures"),
             syntax_case("aromatic lowercase naphthalene"),
             compound_smiles("Sodium Acetate"),
             syntax_case("quadruple bond symbol"),
             syntax_case("aromatic bond between aromatic atoms"),
             syntax_case("wildcard atom"),
         }) {
        text += smiles + "\n";
    }
    EXPECT_EQ(parse_ok(text).size(), 10U);
}

TEST(CandidateList, SharesTheLineRules) {
    const auto ethanol = compound_smiles("Ethanol");
    const auto ethane = compound_smiles("Ethane");
    const auto candidates =
        parse_ok("\xEF\xBB\xBF# list\r\n\r\n  " + ethanol + " ethanol\r\n\r\n" + ethane + "\r\n");
    const std::vector<Candidate> expected{
        {.smiles = ethanol, .name = "ethanol", .line = 3},
        {.smiles = ethane, .name = "", .line = 5},
    };
    EXPECT_EQ(candidates, expected);
}

TEST(CandidateList, KeepsDuplicates) {
    // Removing duplicates needs canonical SMILES, which is chem's job.
    const auto ethanol = compound_smiles("Ethanol");
    EXPECT_EQ(parse_ok(ethanol + " a\n" + ethanol + " a\n").size(), 2U);
}

// ---- Errors -------------------------------------------------------------------

TEST(CandidateList, EmptyListIsError) {
    expect_error("", {"no candidates"});
    expect_error("# only a comment\n\n", {"no candidates"});
}

TEST(CandidateList, CommaSuggestsACsvFile) {
    const auto ethanol = compound_smiles("Ethanol");
    expect_error(ethanol + " ethanol\n" + ethanol + ",ethanol\n", {"line 2", "','", "CSV"});
}

TEST(CandidateList, CharactersOutsideTheSmilesAlphabetAreErrors) {
    // Ethanol's SMILES with one character outside the SMILES alphabet inserted.
    const auto ethanol = compound_smiles("Ethanol");
    expect_error("\"" + ethanol + "\" ethanol\n", {"line 1", "'\"'"});
    expect_error(with_inserted(ethanol, 2, ";") + "\n", {"line 1", "';'"});
    expect_error(with_inserted(ethanol, 2, "\xC3\xA4") + "\n", {"line 1", "non-ASCII"});
    // 0x80 is the first non-ASCII byte; the message must not echo the raw byte.
    expect_error(with_inserted(ethanol, 1, "\x80") + "\n", {"line 1", "non-ASCII"});
    expect_error(with_inserted(ethanol, 1, "\x01") + "\n", {"line 1", "control"});
}

TEST(CandidateList, HeaderLineIsError) {
    const auto ethanol = compound_smiles("Ethanol");
    expect_error("SMILES Name\n" + ethanol + " ethanol\n", {"line 1", "header"});
    expect_error("smiles\n" + ethanol + "\n", {"line 1", "header"});
}

}  // namespace
