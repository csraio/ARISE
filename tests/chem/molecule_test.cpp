#include "arise/chem/molecule.hpp"

#include "arise/core/error.hpp"
#include "arise/core/formula.hpp"
#include "arise/io/pubchem.hpp"
#include "pubchem_fixture.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::test::compound_smiles;
using arise::test::syntax_case;

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
    const auto molecule = parse_ok(compound_smiles("Ethyl Acetate"));
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
    const auto molecule = parse_ok(compound_smiles("Benzene"));
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{1, 1, 1, 1, 1, 1}));
    ASSERT_EQ(molecule.bonds().size(), 6U);
    int doubles = 0;
    for (const auto& bond : molecule.bonds()) {
        doubles += bond.order == BondOrder::Double ? 1 : 0;
    }
    EXPECT_EQ(doubles, 3);
}

TEST(Smiles, ReadsTripleBonds) {
    const auto molecule = parse_ok(compound_smiles("Acetonitrile"));  // acetonitrile
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{3, 0, 0}));
    EXPECT_EQ(molecule.bonds()[1].order, BondOrder::Triple);
}

TEST(Smiles, ReadsNestedBranches) {
    const auto molecule = parse_ok(compound_smiles("Neopentane"));  // neopentane
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{3, 0, 3, 3, 3}));
    EXPECT_EQ(molecule.bonds_of(1).size(), 4U);
}

TEST(Smiles, ReadsTwoDigitRingClosures) {
    const auto molecule = parse_ok(syntax_case("two-digit ring closures"));  // cyclohexane
    EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{2, 2, 2, 2, 2, 2}));
    EXPECT_EQ(molecule.bonds().size(), 6U);
}

TEST(Smiles, TwoDigitRingNumbersAreDistinct) {
    // %12 and %21 are different rings: 12 joins atoms 0-4, 21 joins atoms 2-6.
    const auto molecule = parse_ok(syntax_case("two-digit ring numbers 12 and 21"));
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
    for (const auto& smiles : {
             syntax_case("ring bond order on the opening digit"),
             syntax_case("ring bond order on the closing digit"),
             syntax_case("ring bond order on both digits"),
         }) {
        const auto molecule = parse_ok(smiles);  // cyclohexene
        EXPECT_EQ(hydrogen_counts(molecule), (std::vector<int>{1, 2, 2, 2, 2, 1})) << smiles;
    }
}

TEST(Smiles, ImplicitHydrogensFollowNormalValences) {
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Methylamine"))),
              (std::vector<int>{3, 2}));  // N 3
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Formaldehyde"))),
              (std::vector<int>{2, 0}));                                       // O 2
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Dimethyl Sulfone"))),  // S 6
              (std::vector<int>{3, 0, 0, 0, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Trimethylphosphine"))),
              (std::vector<int>{3, 0, 3, 3}));  // P 3
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Dichloromethane"))),
              (std::vector<int>{0, 2, 0}));
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Bromofluoroiodomethane"))),
              (std::vector<int>{0, 1, 0, 0}));
    EXPECT_EQ(hydrogen_counts(parse_ok("B")), (std::vector<int>{3}));  // borane
}

TEST(Smiles, ReadsBracketAtoms) {
    const auto labelled = parse_ok(compound_smiles("Carbon-13C"));
    ASSERT_EQ(labelled.atoms().size(), 1U);
    EXPECT_EQ(
        labelled.atoms()[0],
        (Atom{.element = Element::C, .charge = 0, .hydrogens = 4, .isotope = 13, .bracket = true}));

    const auto silane = parse_ok(compound_smiles("Tetramethylsilane"));  // tetramethylsilane
    EXPECT_EQ(silane.atoms()[1].element, Element::Si);
    EXPECT_EQ(hydrogen_counts(silane), (std::vector<int>{3, 0, 3, 3, 3}));

    const auto ammonium = parse_ok(compound_smiles("Ammonium"));
    EXPECT_EQ(ammonium.atoms()[0].charge, 1);
    EXPECT_EQ(ammonium.atoms()[0].hydrogens, 4);
    EXPECT_EQ(parse_ok(compound_smiles("Methoxide")).atoms()[1].charge, -1);
    EXPECT_EQ(parse_ok(syntax_case("doubly charged CH2 written ++")).atoms()[0].charge, 2);
    EXPECT_EQ(parse_ok(syntax_case("doubly charged CH2 written +2")).atoms()[0].charge, 2);
}

TEST(Smiles, BracketAtomsGetNoImplicitHydrogens) {
    // [CH] has exactly one hydrogen, whatever its bonds; [C] has none.
    EXPECT_EQ(hydrogen_counts(parse_ok(syntax_case("CH in brackets between methyls"))),
              (std::vector<int>{3, 1, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok(syntax_case("bare carbon in brackets"))),
              (std::vector<int>{0}));
}

TEST(Smiles, ExplicitHydrogenAtomsCountTowardTheirNeighbour) {
    const auto methane = parse_ok(syntax_case("explicit hydrogen atoms on methane"));
    ASSERT_EQ(methane.atoms().size(), 5U);
    EXPECT_EQ(methane.atoms()[1].hydrogens, 0);  // four bonds already: no implicit H
    EXPECT_EQ(methane.hydrogen_count(1), 4);
    EXPECT_EQ(methane.valence(1), 4);
}

TEST(Smiles, StereoMarksAreAcceptedAndIgnored) {
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("(-)-2-Butanol"))),
              (std::vector<int>{3, 1, 1, 2, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("(+)-2-Butanol"))),
              (std::vector<int>{3, 1, 1, 2, 3}));
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("trans-1,2-Difluoroethylene"))),
              (std::vector<int>{0, 1, 1, 0}));
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("cis-1,2-Difluoroethylene"))),
              (std::vector<int>{0, 1, 1, 0}));
}

TEST(Smiles, NeighboursAndValence) {
    const auto molecule = parse_ok(compound_smiles("Acetic Acid"));  // acetic acid
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
    expect_error(syntax_case("aromatic lowercase atom after a Kekule ring"), {"aromatic"});
    expect_error(syntax_case("aromatic bond symbol"), {"column 3", "aromatic"});
}

TEST(Smiles, UnsupportedElementsAndNotationAreErrors) {
    expect_error(syntax_case("sodium"), {"column 3", "Na"});
    expect_error("CSi", {"column 2", "bracket"});  // Si needs brackets: [Si]
    expect_error(syntax_case("quadruple bond symbol"), {"column 2", "'$'"});
    expect_error(syntax_case("atom class"), {"class"});
    expect_error(syntax_case("tetrahedral stereo class"), {"column 4", "stereo"});
}

TEST(Smiles, MoreThanOneComponentIsError) {
    expect_error("CC.O", {"column 3", "component"});
}

TEST(Smiles, UnbalancedBranchesAreErrors) {
    expect_error(syntax_case("unclosed branch"), {"branch"});
    expect_error(syntax_case("branch closed without opening"), {"column 3", "branch"});
    expect_error(syntax_case("empty branch"), {"column 3", "empty branch"});
    expect_error(syntax_case("branch before any atom"), {"column 1", "branch"});
}

TEST(Smiles, BadRingClosuresAreErrors) {
    expect_error("C1CC", {"ring 1", "not closed"});
    expect_error("C11", {"column 3", "itself"});
    expect_error(syntax_case("ring bond orders in conflict"), {"column 10", "ring 1", "conflict"});
    expect_error(syntax_case("percent with one digit"), {"column 2", "%"});
    // Two closures between the same two atoms would be a double bond in disguise.
    expect_error("C12CC12", {"duplicate bond"});
    expect_error("C1C1", {"duplicate bond"});  // chain bond plus ring bond
    // A ring opened inside a branch and closed after it, onto the branch root.
    expect_error(syntax_case("ring closure onto the branch root"), {"duplicate bond"});
}

TEST(Smiles, DanglingBondsAreErrors) {
    expect_error(syntax_case("bond at the end"), {"bond"});
    expect_error(syntax_case("bond at the start"), {"column 1", "bond"});
    expect_error(syntax_case("bond into an empty branch"), {"column 3", "bond"});
}

TEST(Smiles, MalformedBracketAtomsAreErrors) {
    expect_error("C[CH4", {"column 2", "']'"});
    expect_error(syntax_case("empty bracket"), {"column 2", "empty"});
    expect_error(syntax_case("aromatic bracket atom"), {"aromatic"});
}

// ---- Structure issues ---------------------------------------------------------------

TEST(StructureIssue, NormalMoleculesHaveNone) {
    for (const auto& smiles : {
             compound_smiles("Ethyl Acetate"),
             compound_smiles("Acetonitrile"),
             compound_smiles("Dimethyl Sulfone"),
             compound_smiles("Benzene"),
             compound_smiles("Carbon-13C"),
             compound_smiles("Tetramethylsilane"),
             compound_smiles("Trimethylphosphine oxide"),
         }) {
        EXPECT_EQ(structure_issue(parse_ok(smiles)), std::nullopt) << smiles;
    }
}

TEST(StructureIssue, ChargedSpeciesAreFlaggedFirst) {
    EXPECT_EQ(structure_issue(parse_ok(compound_smiles("Ammonium"))), StructureIssue::Charged);
    EXPECT_EQ(structure_issue(parse_ok(compound_smiles("Butane-1,4-diolate"))),
              StructureIssue::Charged);
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("cationic CH2"))),
              StructureIssue::Charged);  // not radical
}

TEST(StructureIssue, RadicalsAndCarbenesAreFlagged) {
    // From the real PubChem C4H8O2 result: an O-CH fragment with one hydrogen too few.
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("O-CH fragment one hydrogen short"))),
              StructureIssue::Radical);
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("methoxy radical"))), StructureIssue::Radical);
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("CH2 radical"))),
              StructureIssue::Radical);  // carbene
    // S with three bonds: valence 3, between its normal valences 2 and 4.
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("sulfur with three bonds"))),
              StructureIssue::Radical);
    // ...whereas one more hydrogen makes valence 4, which is normal for S.
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("sulfur with three bonds and one hydrogen"))),
              std::nullopt);
}

TEST(StructureIssue, ExcessValenceIsFlagged) {
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("carbon with five bonds"))),
              StructureIssue::UnusualValence);
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("CH5"))), StructureIssue::UnusualValence);
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("CF5"))), StructureIssue::UnusualValence);
}

TEST(StructureIssue, NeutralNitrogenWithFiveSingleBondsIsUnusual) {
    // Valence 5 needs a double bond (N(=O)=O); five sigma bonds cannot exist on neutral N
    // (user decision 2026-11-01; checked against RDKit by Provenance.Valences*).
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("neutral nitrogen with four single bonds"))),
              StructureIssue::UnusualValence);
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("neutral nitrogen with five single bonds"))),
              StructureIssue::UnusualValence);
    EXPECT_EQ(structure_issue(parse_ok(compound_smiles("Nitromethane (uncharged nitro)"))),
              std::nullopt);
}

TEST(StructureIssue, IodineMayBeHypervalent) {
    // Iodine valences 1, 3, 5 (hypervalent iodine reagents; user decision 2026-11-01).
    EXPECT_EQ(structure_issue(parse_ok(compound_smiles("Iodobenzene diacetate"))), std::nullopt);
    EXPECT_EQ(hydrogen_counts(parse_ok(compound_smiles("Trimethyl-lambda3-iodane"))),
              (std::vector<int>{3, 0, 3, 3}));
    EXPECT_EQ(structure_issue(parse_ok(syntax_case("iodine with two bonds in brackets"))),
              StructureIssue::Radical);
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

// The formula's degrees of unsaturation count every C as 4-valent and every O as
// 2-valent. An atom away from that valence (a charge pair in a charge-separated structure,
// such as a CH2+ and an O- in one molecule) shifts it: DoU = rings + pi + sum(lowest valence -
// valence) / 2.
int formula_unsaturation(const Molecule& molecule) {
    int excess = 0;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        const auto element = molecule.atoms()[i].element;
        const int lowest = element == arise::chem::Element::C ? 4 : 2;  // C4H8O2: C, O
        excess += lowest - molecule.valence(i);
    }
    return graph_unsaturation(molecule) + (excess / 2);
}

// One PubChem candidate against the formula it was searched with.
void expect_matches_formula(const arise::io::Candidate& candidate,
                            const std::map<std::string, int>& expected_counts,
                            int expected_unsaturation) {
    const auto molecule = parse_smiles(candidate.smiles);
    ASSERT_TRUE(molecule.ok()) << candidate.smiles;
    EXPECT_EQ(element_counts(molecule.value()), expected_counts) << candidate.smiles;
    if (!structure_issue(molecule.value()).has_value()) {
        // Valid closed-shell structures: rings + pi bonds, corrected for atoms away from
        // their lowest valence, equal the formula's DoU.
        EXPECT_EQ(formula_unsaturation(molecule.value()), expected_unsaturation)
            << candidate.smiles;
    }
}

TEST(Smiles, EveryPubChemC4H8O2StructureMatchesTheSearchedFormula) {
    // PubChem returned these for the formula C4H8O2, so every parsed molecule must
    // contain exactly 4 C, 8 H and 2 O. The reference comes from PubChem, not from
    // ARISE's own valence rules, so this checks the reader on all 144 real molecules.
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
    // After io cleaning: 144 structures. Counted by RDKit (radical electrons): exactly 9
    // are radical/carbene fragments ([CH], [CH2]); the other 135 are closed-shell
    // candidates, 13 of them charge-separated.
    const auto imported = arise::io::import_pubchem_csv(arise::test::kPubChemC4H8O2);
    ASSERT_TRUE(imported.ok()) << imported.error().message;
    ASSERT_EQ(imported.value().candidates.size(), 144U);
    const auto tally = tally_issues(imported.value().candidates);
    EXPECT_EQ(tally.radicals, 9);
    EXPECT_EQ(tally.clean, 135);
    EXPECT_EQ(tally.other, 0);
}

}  // namespace
