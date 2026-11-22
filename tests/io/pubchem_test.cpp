#include "arise/io/pubchem.hpp"

#include "arise/core/error.hpp"
#include "arise/core/formula.hpp"
#include "pubchem_fixture.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Formula;
using arise::io::import_pubchem_csv;
using arise::io::PubChemImport;
using arise::io::RemovalReason;
using arise::test::compound_cid;
using arise::test::compound_smiles;
using arise::test::syntax_case;

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

// One CSV row as PubChem writes it: CID, then the quoted SMILES.
std::string row(long long cid, std::string_view smiles, std::string_view end = "\n") {
    return std::to_string(cid) + ",\"" + std::string{smiles} + "\"" + std::string{end};
}

// ---- The search address --------------------------------------------------------

TEST(PubChem, FormulaSearchUrlUsesHillNotation) {
    const auto formula = Formula::parse("O2H8C4");
    ASSERT_TRUE(formula.ok());
    EXPECT_EQ(arise::io::pubchem_formula_search_url(formula.value()),
              "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/C4H8O2/"
              "property/ConnectivitySMILES/CSV");
}

// ---- A real response: C4H8O2, recorded 2026-10-31 ---------------------------------
// 494 rows; 63 mixtures; 10 ions (charges not summing to zero) of the rest; 421 left,
// 144 distinct. Every record's mixture/ion/kept verdict is checked against RDKit by
// Provenance.PubChemCleaningRemovesExactlyMixturesAndIons.

TEST(PubChem, RecordedC4H8O2ResponseIsCleanedAsExpected) {
    const auto imported = import_ok(arise::test::kPubChemC4H8O2);
    EXPECT_EQ(imported.rows, 494U);
    EXPECT_EQ(imported.removed_count(RemovalReason::Mixture), 63U);
    EXPECT_EQ(imported.removed_count(RemovalReason::Charged), 10U);
    EXPECT_EQ(imported.removed_count(RemovalReason::DuplicateConnectivity), 421U - 144U);
    EXPECT_EQ(imported.candidates.size(), 144U);
    EXPECT_EQ(imported.candidates.size() + imported.removed.size(), imported.rows);
}

TEST(PubChem, RecordedResponseKeepsTheLowestCidOfEachStructure) {
    const auto imported = import_ok(arise::test::kPubChemC4H8O2);
    const auto ethyl_acetate =
        std::find_if(imported.candidates.begin(), imported.candidates.end(),
                     [](const arise::io::Candidate& candidate) {
                         return candidate.smiles == compound_smiles("Ethyl Acetate");
                     });
    ASSERT_NE(ethyl_acetate, imported.candidates.end());
    EXPECT_EQ(ethyl_acetate->name, "CID " + std::to_string(compound_cid("Ethyl Acetate")));
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
    // Listed in descending CID order; kept in ascending order.
    const std::pair ester{compound_cid("Ethyl Acetate"), compound_smiles("Ethyl Acetate")};
    const std::pair acid{compound_cid("Butyric Acid"), compound_smiles("Butyric Acid")};
    const auto& [low, high] = std::minmax(ester, acid);
    const auto imported =
        import_ok(std::string{kHeader} + row(high.first, high.second) + row(low.first, low.second));
    ASSERT_EQ(imported.candidates.size(), 2U);
    EXPECT_EQ(imported.candidates[0].smiles, low.second);
    EXPECT_EQ(imported.candidates[0].name, "CID " + std::to_string(low.first));
    EXPECT_EQ(imported.candidates[0].line, 3U);
    EXPECT_EQ(imported.candidates[1].name, "CID " + std::to_string(high.first));
    EXPECT_EQ(imported.candidates[1].line, 2U);
}

TEST(PubChem, DuplicateKeepsTheLowestCidWhateverTheRowOrder) {
    const auto ethanol = compound_smiles("Ethanol");
    const auto imported =
        import_ok(std::string{kHeader} + row(500, ethanol) + row(7, ethanol) + row(90, ethanol));
    ASSERT_EQ(imported.candidates.size(), 1U);
    EXPECT_EQ(imported.candidates[0].name, "CID 7");
    const std::vector<arise::io::Removal> expected{
        {.cid = 90, .smiles = ethanol, .reason = RemovalReason::DuplicateConnectivity},
        {.cid = 500, .smiles = ethanol, .reason = RemovalReason::DuplicateConnectivity},
    };
    EXPECT_EQ(imported.removed, expected);
}

TEST(PubChem, RemovesMixturesAndIons) {
    const auto imported =
        import_ok(std::string{kHeader} +
                  row(1, syntax_case("mixture with a charged component")) +  // also charged
                  row(2, compound_smiles("Butane-1,4-diolate")) +            // net charge -2
                  row(3, compound_smiles("Methylammonium")) +                // net charge +1
                  row(4, compound_smiles("Ethanol")));
    ASSERT_EQ(imported.candidates.size(), 1U);
    EXPECT_EQ(imported.candidates[0].smiles, compound_smiles("Ethanol"));
    EXPECT_EQ(imported.removed_count(RemovalReason::Mixture), 1U);  // mixture wins
    EXPECT_EQ(imported.removed_count(RemovalReason::Charged), 2U);
}

TEST(PubChem, KeepsChargeSeparatedNeutralMolecules) {
    // PubChem writes nitro groups, N-oxides, zwitterions and isocyanides with formal
    // charges that sum to zero; they are real neutral compounds (user decision 2026-11-01).
    const auto imported = import_ok(std::string{kHeader} + row(1, compound_smiles("Nitroethane")) +
                                    row(2, compound_smiles("Pyridine N-Oxide")) +
                                    row(3, compound_smiles("alpha-Glycine")) +
                                    row(4, compound_smiles("Methyl isocyanide")));
    EXPECT_EQ(imported.candidates.size(), 4U);
    EXPECT_TRUE(imported.removed.empty());
}

TEST(PubChem, ChargesAreSummedInEveryNotation) {
    // '++' and '+2' are both +2; one -1 per '-'.
    const auto imported =
        import_ok(std::string{kHeader} +
                  row(1, syntax_case("CH2 dication and two oxide anions, written ++")) +  // kept
                  row(2, syntax_case("CH2 dication and two oxide anions, written +2")) +  // kept
                  row(3, syntax_case("CH2 dication and one oxide anion")) +               // ion
                  row(4, syntax_case("oxide dianion written --")) +                       // ion
                  row(5, syntax_case("oxide dianion written -2, two cations")));          // kept
    EXPECT_EQ(imported.candidates.size(), 3U);
    EXPECT_EQ(imported.removed_count(RemovalReason::Charged), 2U);
}

TEST(PubChem, ExplicitSingleBondIsNotACharge) {
    // '-' outside brackets is a bond, not a charge.
    const auto imported =
        import_ok(std::string{kHeader} + row(1, syntax_case("explicit single bond")) +
                  row(2, compound_smiles("Carbon-13C")));
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
    const auto imported = import_ok("\"CID\",\"ConnectivitySMILES\"\r\n" +
                                    row(7, compound_smiles("Ethanol"), "\r\n") +
                                    row(8, compound_smiles("Ethane"), ""));
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
    expect_error("\"CID\",\"IsomericSMILES\"\n" + row(7, compound_smiles("Ethanol")),
                 {"line 1", "header"});
}

TEST(PubChem, BadCidIsError) {
    const auto field = ",\"" + compound_smiles("Ethanol") + "\"\n";
    expect_error(std::string{kHeader} + "abc" + field, {"line 2", "CID"});
    expect_error(std::string{kHeader} + "0" + field, {"line 2", "CID"});
    expect_error(std::string{kHeader} + field, {"line 2", "CID"});
}

TEST(PubChem, MalformedSmilesFieldIsError) {
    const auto ethanol = compound_smiles("Ethanol");
    expect_error(std::string{kHeader} + "7," + ethanol + "\n", {"line 2", "quoted"});
    expect_error(std::string{kHeader} + "7,\"" + ethanol + "\n", {"line 2", "quoted"});
    expect_error(std::string{kHeader} + "7\n", {"line 2"});
    expect_error(std::string{kHeader} + "7,\"\"\n", {"line 2", "empty"});
    // Ethanol's SMILES with a character outside the SMILES alphabet inserted.
    expect_error(std::string{kHeader} + row(7, std::string{ethanol}.insert(2, ";")),
                 {"line 2", "';'"});
}

}  // namespace
