#include "arise/io/pubchem.hpp"

#include "arise/core/error.hpp"
#include "arise/core/formula.hpp"
#include "pubchem_fixture.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Formula;
using arise::io::import_pubchem_csv;
using arise::io::PubChemImport;
using arise::io::RemovalReason;

PubChemImport import_ok(std::string_view text) {
    auto result = import_pubchem_csv(text);
    EXPECT_TRUE(result.ok()) << (result.ok() ? "" : result.error().message);
    return result.ok() ? std::move(result).value() : PubChemImport{};
}

void expect_error(std::string_view text, const std::vector<std::string>& fragments) {
    const auto result = import_pubchem_csv(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, ErrorCode::ParseError) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

constexpr std::string_view kHeader = "\"CID\",\"ConnectivitySMILES\"\n";

// ---- The search address --------------------------------------------------------

TEST(PubChem, FormulaSearchUrlUsesHillNotation) {
    const auto formula = Formula::parse("O2H8C4");
    ASSERT_TRUE(formula.ok());
    EXPECT_EQ(arise::io::pubchem_formula_search_url(formula.value()),
              "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/C4H8O2/"
              "property/ConnectivitySMILES/CSV");
}

// ---- A real response: C4H8O2, recorded 2026-10-31 ---------------------------------
// Expected counts were computed independently with grep/sort on the same file:
// 494 rows; 63 mixtures; 23 charged (of the rest); 408 left, 131 distinct.

TEST(PubChem, RecordedC4H8O2ResponseIsCleanedAsExpected) {
    const auto imported = import_ok(arise::test::kPubChemC4H8O2);
    EXPECT_EQ(imported.rows, 494U);
    EXPECT_EQ(imported.removed_count(RemovalReason::Mixture), 63U);
    EXPECT_EQ(imported.removed_count(RemovalReason::Charged), 23U);
    EXPECT_EQ(imported.removed_count(RemovalReason::DuplicateConnectivity), 408U - 131U);
    EXPECT_EQ(imported.candidates.size(), 131U);
    EXPECT_EQ(imported.candidates.size() + imported.removed.size(), imported.rows);
}

TEST(PubChem, RecordedResponseKeepsTheLowestCidOfEachStructure) {
    const auto imported = import_ok(arise::test::kPubChemC4H8O2);
    const auto ethyl_acetate = std::find_if(
        imported.candidates.begin(), imported.candidates.end(),
        [](const arise::io::Candidate& candidate) { return candidate.smiles == "CCOC(=O)C"; });
    ASSERT_NE(ethyl_acetate, imported.candidates.end());
    EXPECT_EQ(ethyl_acetate->name, "CID 8857");
}

TEST(PubChem, RecordedResponseCandidatesAreDistinctAndClean) {
    const auto imported = import_ok(arise::test::kPubChemC4H8O2);
    std::set<std::string> seen;
    for (const auto& candidate : imported.candidates) {
        EXPECT_TRUE(seen.insert(candidate.smiles).second) << "duplicate " << candidate.smiles;
        EXPECT_EQ(candidate.smiles.find('.'), std::string::npos) << candidate.smiles;
    }
}

// ---- Cleaning rules on small hand-made inputs -------------------------------------

TEST(PubChem, KeepsCandidatesInCidOrderWithNamesAndLines) {
    const auto imported =
        import_ok(std::string{kHeader} + "8857,\"CCOC(=O)C\"\n264,\"CCCC(=O)O\"\n");
    ASSERT_EQ(imported.candidates.size(), 2U);
    EXPECT_EQ(imported.candidates[0].smiles, "CCCC(=O)O");
    EXPECT_EQ(imported.candidates[0].name, "CID 264");
    EXPECT_EQ(imported.candidates[0].line, 3U);
    EXPECT_EQ(imported.candidates[1].name, "CID 8857");
    EXPECT_EQ(imported.candidates[1].line, 2U);
}

TEST(PubChem, DuplicateKeepsTheLowestCidWhateverTheRowOrder) {
    const auto imported = import_ok(std::string{kHeader} + "500,\"CCO\"\n7,\"CCO\"\n90,\"CCO\"\n");
    ASSERT_EQ(imported.candidates.size(), 1U);
    EXPECT_EQ(imported.candidates[0].name, "CID 7");
    const std::vector<arise::io::Removal> expected{
        {.cid = 90, .smiles = "CCO", .reason = RemovalReason::DuplicateConnectivity},
        {.cid = 500, .smiles = "CCO", .reason = RemovalReason::DuplicateConnectivity},
    };
    EXPECT_EQ(imported.removed, expected);
}

TEST(PubChem, RemovesMixturesAndChargedSpecies) {
    const auto imported = import_ok(std::string{kHeader} +
                                    "1,\"CC(C)O.[C-]#[O+]\"\n"  // mixture, also charged
                                    "2,\"C(CC[O-])C[O-]\"\n"    // charged
                                    "3,\"C[NH3+]\"\n"           // charged
                                    "4,\"CCO\"\n");
    ASSERT_EQ(imported.candidates.size(), 1U);
    EXPECT_EQ(imported.candidates[0].smiles, "CCO");
    EXPECT_EQ(imported.removed_count(RemovalReason::Mixture), 1U);  // mixture wins
    EXPECT_EQ(imported.removed_count(RemovalReason::Charged), 2U);
}

TEST(PubChem, ExplicitSingleBondIsNotACharge) {
    // '-' outside brackets is a bond, not a charge.
    const auto imported = import_ok(std::string{kHeader} + "1,\"CC-O\"\n2,\"[13CH4]\"\n");
    EXPECT_EQ(imported.candidates.size(), 2U);
    EXPECT_TRUE(imported.removed.empty());
}

TEST(PubChem, HeaderOnlyIsAnEmptyResult) {
    // No compound has this formula: valid, nothing to prune.
    const auto imported = import_ok(kHeader);
    EXPECT_TRUE(imported.candidates.empty());
    EXPECT_EQ(imported.rows, 0U);
}

TEST(PubChem, AcceptsCrLfAndAMissingFinalNewline) {
    const auto imported = import_ok("\"CID\",\"ConnectivitySMILES\"\r\n7,\"CCO\"\r\n8,\"CC\"");
    EXPECT_EQ(imported.candidates.size(), 2U);
}

TEST(PubChem, ReasonsHaveReadableNames) {
    EXPECT_EQ(arise::io::to_string(RemovalReason::Mixture), "mixture");
    EXPECT_EQ(arise::io::to_string(RemovalReason::Charged), "charged");
    EXPECT_EQ(arise::io::to_string(RemovalReason::DuplicateConnectivity), "duplicate structure");
}

// ---- Not a PubChem CSV ------------------------------------------------------------

TEST(PubChem, EmptyTextIsError) {
    expect_error("", {"header"});
}

TEST(PubChem, WrongHeaderIsError) {
    // e.g. an HTML error page or a different property list.
    expect_error("<html><body>Service unavailable</body></html>\n", {"line 1", "header"});
    expect_error("\"CID\",\"IsomericSMILES\"\n7,\"CCO\"\n", {"line 1", "header"});
}

TEST(PubChem, BadCidIsError) {
    expect_error(std::string{kHeader} + "abc,\"CCO\"\n", {"line 2", "CID"});
    expect_error(std::string{kHeader} + "0,\"CCO\"\n", {"line 2", "CID"});
    expect_error(std::string{kHeader} + ",\"CCO\"\n", {"line 2", "CID"});
}

TEST(PubChem, MalformedSmilesFieldIsError) {
    expect_error(std::string{kHeader} + "7,CCO\n", {"line 2", "quoted"});
    expect_error(std::string{kHeader} + "7,\"CCO\n", {"line 2", "quoted"});
    expect_error(std::string{kHeader} + "7\n", {"line 2"});
    expect_error(std::string{kHeader} + "7,\"\"\n", {"line 2", "empty"});
    expect_error(std::string{kHeader} + "7,\"CC;O\"\n", {"line 2", "';'"});
}

}  // namespace
