#include "arise/chem/molecule.hpp"

#include "arise/core/error.hpp"
#include "arise/core/formula.hpp"
#include "arise/io/pubchem.hpp"
#include "pubchem_fixture.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::chem::Atom;
using arise::chem::Bond;
using arise::chem::BondOrder;
using arise::chem::Element;
using arise::chem::Molecule;
using arise::chem::parse_smiles;
using arise::chem::structure_issue;
using arise::chem::StructureIssue;
using arise::core::ErrorCode;

Molecule parse_ok(std::string_view smiles) {
    auto result = parse_smiles(smiles);
    EXPECT_TRUE(result.ok()) << smiles << ": " << (result.ok() ? "" : result.error().message);
    return result.ok() ? std::move(result).value() : Molecule{};
}

void expect_error(std::string_view smiles, const std::vector<std::string>& fragments) {
    const auto result = parse_smiles(smiles);
    ASSERT_FALSE(result.ok()) << "accepted '" << smiles << "'";
    EXPECT_EQ(result.error().code, ErrorCode::ParseError) << result.error().message;
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

std::vector<int> hydrogen_counts(const Molecule& molecule) {
    std::vector<int> counts;
    counts.reserve(molecule.atoms().size());
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        counts.push_back(molecule.hydrogen_count(i));
    }
    return counts;
}

std::vector<Element> elements(const Molecule& molecule) {
    std::vector<Element> result;
    result.reserve(molecule.atoms().size());
    for (const auto& atom : molecule.atoms()) {
        result.push_back(atom.element);
    }
    return result;
}

// ---- Reading molecules ---------------------------------------------------------

TEST(Smiles, ReadsEthylAcetate) {
    // CH3-CH2-O-C(=O)-CH3, as PubChem writes it.
    const auto molecule = parse_ok("CCOC(=O)C");
    const std::vector<Element> expected_elements{
        Element::C, Element::C, Element::O, Element::C, Element::O, Element::C,
    };
    EXPECT_EQ(elements(molecule), expected_elements);
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{3, 2, 0, 0, 0, 3}));
    const std::vector<Bond> expected_bonds{
        {.first = 0, .second = 1, .order = BondOrder::Single},
        {.first = 1, .second = 2, .order = BondOrder::Single},
        {.first = 2, .second = 3, .order = BondOrder::Single},
        {.first = 3, .second = 4, .order = BondOrder::Double},
        {.first = 3, .second = 5, .order = BondOrder::Single},
    };
    EXPECT_EQ(std::vector<Bond>(molecule.bonds().begin(), molecule.bonds().end()), expected_bonds);
}

TEST(Smiles, ReadsKekuleBenzeneAsSixMethines) {
    const auto molecule = parse_ok("C1=CC=CC=C1");
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{1, 1, 1, 1, 1, 1}));
    ASSERT_EQ(molecule.bonds().size(), 6U);
    int doubles = 0;
    for (const auto& bond : molecule.bonds()) {
        doubles += bond.order == BondOrder::Double ? 1 : 0;
    }
    EXPECT_EQ(doubles, 3);
}

TEST(Smiles, ReadsTripleBonds) {
    const auto molecule = parse_ok("CC#N");  // acetonitrile
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{3, 0, 0}));
    EXPECT_EQ(molecule.bonds()[1].order, BondOrder::Triple);
}

TEST(Smiles, ReadsNestedBranches) {
    const auto molecule = parse_ok("CC(C)(C)C");  // neopentane
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{3, 0, 3, 3, 3}));
    EXPECT_EQ(molecule.bonds_of(1).size(), 4U);
}

TEST(Smiles, ReadsTwoDigitRingClosures) {
    const auto molecule = parse_ok("C%10CCCCC%10");  // cyclohexane
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{2, 2, 2, 2, 2, 2}));
    EXPECT_EQ(molecule.bonds().size(), 6U);
}

TEST(Smiles, TwoDigitRingNumbersAreDistinct) {
    // %12 and %21 are different rings: 12 joins atoms 0-4, 21 joins atoms 2-6.
    const auto molecule = parse_ok("C%12CC%21CC%12CC%21");
    std::vector<std::pair<std::size_t, std::size_t>> ring_bonds;
    for (const auto& bond : molecule.bonds()) {
        if (bond.second != bond.first + 1) {
            ring_bonds.emplace_back(bond.first, bond.second);
        }
    }
    const std::vector<std::pair<std::size_t, std::size_t>> expected{{0, 4}, {2, 6}};
    EXPECT_EQ(ring_bonds, expected);
}

TEST(Smiles, RingBondOrderMayBeWrittenOnEitherEnd) {
    for (const auto* smiles : {"C=1CCCCC1", "C1CCCCC=1", "C=1CCCCC=1"}) {
        const auto molecule = parse_ok(smiles);  // cyclohexene
        EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{1, 2, 2, 2, 2, 1})) << smiles;
    }
}

TEST(Smiles, ImplicitHydrogensFollowNormalValences) {
    EXPECT_EQ(hydrogen_counts(parse_ok("CN")), (std::vector<int>{3, 2}));   // N 3
    EXPECT_EQ(hydrogen_counts(parse_ok("C=O")), (std::vector<int>{2, 0}));  // O 2
    EXPECT_EQ(hydrogen_counts(parse_ok("CS(=O)(=O)C")),                     // S 6
              (std::vector<int>{3, 0, 0, 0, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok("CP(C)C")), (std::vector<int>{3, 0, 3, 3}));  // P 3
    EXPECT_EQ(hydrogen_counts(parse_ok("ClCCl")), (std::vector<int>{0, 2, 0}));
    EXPECT_EQ(hydrogen_counts(parse_ok("BrC(F)I")), (std::vector<int>{0, 1, 0, 0}));
    EXPECT_EQ(hydrogen_counts(parse_ok("B")), (std::vector<int>{3}));  // borane
}

TEST(Smiles, ReadsBracketAtoms) {
    const auto labelled = parse_ok("[13CH4]");
    ASSERT_EQ(labelled.atoms().size(), 1U);
    EXPECT_EQ(
        labelled.atoms()[0],
        (Atom{.element = Element::C, .charge = 0, .hydrogens = 4, .isotope = 13, .bracket = true}));

    const auto silane = parse_ok("C[Si](C)(C)C");  // tetramethylsilane
    EXPECT_EQ(silane.atoms()[1].element, Element::Si);
    EXPECT_EQ(hydrogen_counts(silane), (std::vector<int>{3, 0, 3, 3, 3}));

    const auto ammonium = parse_ok("[NH4+]");
    EXPECT_EQ(ammonium.atoms()[0].charge, 1);
    EXPECT_EQ(ammonium.atoms()[0].hydrogens, 4);
    EXPECT_EQ(parse_ok("C[O-]").atoms()[1].charge, -1);
    EXPECT_EQ(parse_ok("[CH2++]").atoms()[0].charge, 2);
    EXPECT_EQ(parse_ok("[CH2+2]").atoms()[0].charge, 2);
}

TEST(Smiles, BracketAtomsGetNoImplicitHydrogens) {
    // [CH] has exactly one hydrogen, whatever its bonds; [C] has none.
    EXPECT_EQ(hydrogen_counts(parse_ok("C[CH]C")), (std::vector<int>{3, 1, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok("[C]")), (std::vector<int>{0}));
}

TEST(Smiles, ExplicitHydrogenAtomsCountTowardTheirNeighbour) {
    const auto methane = parse_ok("[H]C([H])([H])[H]");
    ASSERT_EQ(methane.atoms().size(), 5U);
    EXPECT_EQ(methane.atoms()[1].hydrogens, 0);  // four bonds already: no implicit H
    EXPECT_EQ(methane.hydrogen_count(1), 4);
    EXPECT_EQ(methane.valence(1), 4);
}

TEST(Smiles, StereoMarksAreAcceptedAndIgnored) {
    EXPECT_EQ(hydrogen_counts(parse_ok("C[C@@H](O)CC")), (std::vector<int>{3, 1, 1, 2, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok("C[C@H](O)CC")), (std::vector<int>{3, 1, 1, 2, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok("F/C=C/F")), (std::vector<int>{0, 1, 1, 0}));
    EXPECT_EQ(hydrogen_counts(parse_ok("F/C=C\\F")), (std::vector<int>{0, 1, 1, 0}));
}

TEST(Smiles, NeighboursAndValence) {
    const auto molecule = parse_ok("CC(=O)O");  // acetic acid
    std::vector<std::size_t> neighbours;
    for (const auto bond : molecule.bonds_of(1)) {
        neighbours.push_back(molecule.neighbor(bond, 1));
    }
    EXPECT_EQ(neighbours, (std::vector<std::size_t>{0, 2, 3}));
    EXPECT_EQ(molecule.valence(1), 4);
    EXPECT_EQ(molecule.valence(2), 2);
}

// ---- Syntax errors ---------------------------------------------------------------

TEST(Smiles, EmptyInputIsError) {
    expect_error("", {"empty"});
}

TEST(Smiles, AromaticLowercaseIsNotSupportedYet) {
    expect_error("c1ccccc1", {"column 1", "aromatic", "Kekul"});
    expect_error("C1=CC=CC=C1[nH]", {"aromatic"});
    expect_error("CC:C", {"column 3", "aromatic"});
}

TEST(Smiles, UnsupportedElementsAndNotationAreErrors) {
    expect_error("C[Na]", {"column 3", "Na"});
    expect_error("CSi", {"column 2", "bracket"});  // Si needs brackets: [Si]
    expect_error("C$C", {"column 2", "'$'"});
    expect_error("[CH4:1]", {"class"});
    expect_error("C[C@TH1H](O)C", {"column 4", "stereo"});
}

TEST(Smiles, MoreThanOneComponentIsError) {
    expect_error("CC.O", {"column 3", "component"});
}

TEST(Smiles, UnbalancedBranchesAreErrors) {
    expect_error("CC(C", {"branch"});
    expect_error("CC)C", {"column 3", "branch"});
    expect_error("C()C", {"column 3", "empty branch"});
    expect_error("(C)C", {"column 1", "branch"});
}

TEST(Smiles, BadRingClosuresAreErrors) {
    expect_error("C1CC", {"ring 1", "not closed"});
    expect_error("C11", {"column 3", "itself"});
    expect_error("C=1CCCCC#1", {"column 10", "ring 1", "conflict"});
    expect_error("C%1C", {"column 2", "%"});
    // Two closures between the same two atoms would be a double bond in disguise.
    expect_error("C12CC12", {"duplicate bond"});
    expect_error("C1C1", {"duplicate bond"});  // chain bond plus ring bond
    // A ring opened inside a branch and closed after it, onto the branch root.
    expect_error("CC(C1C)1", {"duplicate bond"});
}

TEST(Smiles, DanglingBondsAreErrors) {
    expect_error("CC=", {"bond"});
    expect_error("=CC", {"column 1", "bond"});
    expect_error("C(=)C", {"column 3", "bond"});
}

TEST(Smiles, MalformedBracketAtomsAreErrors) {
    expect_error("C[CH4", {"column 2", "']'"});
    expect_error("C[]C", {"column 2", "empty"});
    expect_error("[cH4]", {"aromatic"});
}

// ---- Structure issues ---------------------------------------------------------------

TEST(StructureIssue, NormalMoleculesHaveNone) {
    for (const auto* smiles : {
             "CCOC(=O)C",
             "CC#N",
             "CS(=O)(=O)C",
             "C1=CC=CC=C1",
             "[13CH4]",
             "C[Si](C)(C)C",
             "CP(=O)(C)C",
         }) {
        EXPECT_EQ(structure_issue(parse_ok(smiles)), std::nullopt) << smiles;
    }
}

TEST(StructureIssue, ChargedSpeciesAreFlaggedFirst) {
    EXPECT_EQ(structure_issue(parse_ok("[NH4+]")), StructureIssue::Charged);
    EXPECT_EQ(structure_issue(parse_ok("C(CC[O-])C[O-]")), StructureIssue::Charged);
    EXPECT_EQ(structure_issue(parse_ok("[CH2+]")), StructureIssue::Charged);  // not radical
}

TEST(StructureIssue, RadicalsAndCarbenesAreFlagged) {
    // From the real PubChem C4H8O2 result: an O-CH fragment with one hydrogen too few.
    EXPECT_EQ(structure_issue(parse_ok("CC(CO)O[CH]")), StructureIssue::Radical);
    EXPECT_EQ(structure_issue(parse_ok("C[O]")), StructureIssue::Radical);
    EXPECT_EQ(structure_issue(parse_ok("[CH2]")), StructureIssue::Radical);  // carbene
    // S with three bonds: valence 3, between its normal valences 2 and 4.
    EXPECT_EQ(structure_issue(parse_ok("C[S](C)C")), StructureIssue::Radical);
    // ...whereas one more hydrogen makes valence 4, which is normal for S.
    EXPECT_EQ(structure_issue(parse_ok("C[SH](C)C")), std::nullopt);
}

TEST(StructureIssue, ExcessValenceIsFlagged) {
    EXPECT_EQ(structure_issue(parse_ok("CC(C)(C)(C)C")), StructureIssue::UnusualValence);
    EXPECT_EQ(structure_issue(parse_ok("[CH5]")), StructureIssue::UnusualValence);
    EXPECT_EQ(structure_issue(parse_ok("FC(F)(F)(F)F")), StructureIssue::UnusualValence);
}

TEST(StructureIssue, IssuesHaveReadableNames) {
    EXPECT_EQ(arise::chem::to_string(StructureIssue::Charged), "charged");
    EXPECT_EQ(arise::chem::to_string(StructureIssue::Radical), "radical");
    EXPECT_EQ(arise::chem::to_string(StructureIssue::UnusualValence), "unusual valence");
}

// ---- The real PubChem result ----------------------------------------------------------

// Tally of structure issues over a candidate list; parse failures are test failures.
struct IssueTally {
    int clean = 0;
    int radicals = 0;
    int other = 0;
};

IssueTally tally_issues(const std::vector<arise::io::Candidate>& candidates) {
    IssueTally tally;
    for (const auto& candidate : candidates) {
        const auto molecule = parse_smiles(candidate.smiles);
        if (!molecule.ok()) {
            ADD_FAILURE() << candidate.smiles << ": " << molecule.error().message;
            continue;
        }
        const auto issue = structure_issue(molecule.value());
        if (!issue.has_value()) {
            ++tally.clean;
        } else if (*issue == StructureIssue::Radical) {
            ++tally.radicals;
        } else {
            ++tally.other;
        }
    }
    return tally;
}

// Atom counts of a parsed molecule, hydrogens included (implicit, bracket and [H] atoms).
std::map<std::string, int> element_counts(const Molecule& molecule) {
    std::map<std::string, int> counts;
    for (const auto& atom : molecule.atoms()) {
        counts[std::string{arise::chem::symbol(atom.element)}] += 1;
        counts["H"] += atom.hydrogens;
    }
    return counts;
}

// Rings plus pi bonds, from the graph alone: (bonds - atoms + 1) + sum(order - 1).
// Explicit [H] atoms count as atoms; implicit hydrogens add a bond and an atom each,
// which cancel, so they can be ignored.
int graph_unsaturation(const Molecule& molecule) {
    const auto bonds = static_cast<int>(molecule.bonds().size());
    const auto atoms = static_cast<int>(molecule.atoms().size());
    int pi = 0;
    for (const auto& bond : molecule.bonds()) {
        pi += static_cast<int>(bond.order) - 1;
    }
    return (bonds - atoms + 1) + pi;
}

// One PubChem candidate against the formula it was searched with.
void expect_matches_formula(const arise::io::Candidate& candidate,
                            const std::map<std::string, int>& expected_counts,
                            int expected_unsaturation) {
    const auto molecule = parse_smiles(candidate.smiles);
    ASSERT_TRUE(molecule.ok()) << candidate.smiles;
    EXPECT_EQ(element_counts(molecule.value()), expected_counts) << candidate.smiles;
    if (!structure_issue(molecule.value()).has_value()) {
        // Valid closed-shell structures: rings + pi bonds equal the formula's DoU.
        EXPECT_EQ(graph_unsaturation(molecule.value()), expected_unsaturation) << candidate.smiles;
    }
}

TEST(Smiles, EveryPubChemC4H8O2StructureMatchesTheSearchedFormula) {
    // PubChem returned these for the formula C4H8O2, so every parsed molecule must
    // contain exactly 4 C, 8 H and 2 O. The reference comes from PubChem, not from
    // ARISE's own valence rules, so this checks the reader on all 131 real molecules.
    const auto formula = arise::core::Formula::parse("C4H8O2");
    ASSERT_TRUE(formula.ok());
    const auto unsaturation = formula.value().degrees_of_unsaturation();
    ASSERT_TRUE(unsaturation.ok());
    const std::map<std::string, int> expected{{"C", 4}, {"H", 8}, {"O", 2}};
    const auto imported = arise::io::import_pubchem_csv(arise::test::kPubChemC4H8O2);
    ASSERT_TRUE(imported.ok()) << imported.error().message;
    for (const auto& candidate : imported.value().candidates) {
        expect_matches_formula(candidate, expected, unsaturation.value());
    }
}

TEST(Smiles, ReadsEveryCleanedPubChemC4H8O2Structure) {
    // After io cleaning: 131 structures. Counted independently with grep: exactly 9 are
    // radical/carbene fragments ([CH], [CH2]); the other 122 are real candidates.
    const auto imported = arise::io::import_pubchem_csv(arise::test::kPubChemC4H8O2);
    ASSERT_TRUE(imported.ok()) << imported.error().message;
    ASSERT_EQ(imported.value().candidates.size(), 131U);
    const auto tally = tally_issues(imported.value().candidates);
    EXPECT_EQ(tally.radicals, 9);
    EXPECT_EQ(tally.clean, 122);
    EXPECT_EQ(tally.other, 0);
}

}  // namespace
