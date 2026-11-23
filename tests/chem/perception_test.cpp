#include "arise/chem/perception.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/io/pubchem.hpp"
#include "arise/nmr1d/dept.hpp"
#include "pubchem_fixture.hpp"
#include "test_molecules.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <string_view>
#include <vector>

namespace {

using arise::test::compound_hybridizations;
using arise::test::compound_smiles;
using arise::test::syntax_case;

using arise::chem::aromatic_rings;
using arise::chem::carbon_hybridization;
using arise::chem::carbon_type;
using arise::chem::count_group;
using arise::chem::Element;
using arise::chem::Group;
using arise::chem::Hybridization;
using arise::chem::Molecule;
using arise::chem::symmetry_classes;
using arise::nmr1d::CarbonType;

Molecule parse_ok(std::string_view smiles) {
    auto result = arise::chem::parse_smiles(smiles);
    EXPECT_TRUE(result.ok()) << smiles << ": " << (result.ok() ? "" : result.error().message);
    return result.ok() ? std::move(result).value() : Molecule{};
}

// Sizes of the symmetry classes that contain carbons, sorted: the 13C peak pattern a
// molecule should give if no peaks overlap by accident.
std::vector<int> carbon_class_sizes(std::string_view smiles) {
    const auto molecule = parse_ok(smiles);
    const auto classes = symmetry_classes(molecule);
    std::map<std::size_t, int> sizes;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (molecule.atoms()[i].element == Element::C) {
            ++sizes[classes[i]];
        }
    }
    std::vector<int> result;
    result.reserve(sizes.size());
    for (const auto& [cls, size] : sizes) {
        result.push_back(size);
    }
    std::sort(result.begin(), result.end());
    return result;
}

int groups(std::string_view smiles, Group group) {
    return count_group(parse_ok(smiles), group);
}

// ---- Carbon types --------------------------------------------------------------

TEST(CarbonType, FollowsAttachedHydrogens) {
    const auto molecule = parse_ok(compound_smiles("Ethyl Acetate"));  // ethyl acetate
    EXPECT_EQ(carbon_type(molecule, 0), CarbonType::CH3);
    EXPECT_EQ(carbon_type(molecule, 1), CarbonType::CH2);
    EXPECT_EQ(carbon_type(molecule, 2), std::nullopt);  // oxygen
    EXPECT_EQ(carbon_type(molecule, 3), CarbonType::C);
    EXPECT_EQ(carbon_type(molecule, 5), CarbonType::CH3);
    EXPECT_EQ(carbon_type(parse_ok(compound_smiles("Benzene")), 0), CarbonType::CH);
}

TEST(CarbonType, MethaneHasNoDeptType) {
    EXPECT_EQ(carbon_type(parse_ok("C"), 0), std::nullopt);
}

TEST(CarbonType, ExplicitHydrogenAtomsCount) {
    EXPECT_EQ(carbon_type(parse_ok(syntax_case("explicit hydrogen atoms on propane C1")), 1),
              CarbonType::CH2);
}

// ---- Carbon hybridization --------------------------------------------------------
// Expected values are RDKit's hybridizations of the PubChem compounds
// (test_compounds.tsv, compound_hybridizations), written as 1 sp, 2 sp2, 3 sp3, 0 none.

int hybridization_code(const Molecule& molecule, std::size_t atom) {
    const auto hybridization = carbon_hybridization(molecule, atom);
    if (!hybridization) {
        return 0;
    }
    switch (*hybridization) {
        case Hybridization::Sp:
            return 1;
        case Hybridization::Sp2:
            return 2;
        case Hybridization::Sp3:
            return 3;
    }
    return 0;
}

std::vector<int> hybridization_codes(std::string_view name) {
    const auto molecule = parse_ok(compound_smiles(name));
    std::vector<int> codes;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        codes.push_back(hybridization_code(molecule, i));
    }
    return codes;
}

TEST(CarbonHybridization, NamesEachHybridization) {
    EXPECT_EQ(to_string(Hybridization::Sp), "sp");
    EXPECT_EQ(to_string(Hybridization::Sp2), "sp2");
    EXPECT_EQ(to_string(Hybridization::Sp3), "sp3");
}

TEST(CarbonHybridization, TripleBondsGiveSp) {
    for (const auto* name : {"Acetonitrile", "Hydrogen Cyanide"}) {
        EXPECT_EQ(hybridization_codes(name), compound_hybridizations(name)) << name;
    }
}

TEST(CarbonHybridization, TwoDoubleBondsGiveSp) {
    // The central carbon of an allene and of a ketene.
    for (const auto* name :
         {"2,3-Pentadiene, 2-methyl-", "Methylketene", "Methyl 5-methylhexa-2,3-dienoate"}) {
        const auto expected = compound_hybridizations(name);
        EXPECT_NE(std::find(expected.begin(), expected.end(), 1), expected.end()) << name;
        EXPECT_EQ(hybridization_codes(name), expected) << name;
    }
}

TEST(CarbonHybridization, OneDoubleBondGivesSp2) {
    // Includes the C=N carbons of an oxime and of a diazo compound.
    for (const auto* name : {
             "Acetone",
             "Formaldehyde",
             "Propylene",
             "Acetone, oxime",
             "3-Diazo-2,2,4,4-tetramethylpentane",
         }) {
        EXPECT_EQ(hybridization_codes(name), compound_hybridizations(name)) << name;
    }
}

TEST(CarbonHybridization, KekuleAromaticRingsGiveSp2) {
    for (const auto* name : {"Benzene", "Pyridine", "Furan", "Naphthalene", "Indole"}) {
        EXPECT_EQ(hybridization_codes(name), compound_hybridizations(name)) << name;
    }
}

TEST(CarbonHybridization, SingleBondsOnlyGiveSp3) {
    for (const auto* name : {"Methane", "Propane", "Cyclohexane", "Neopentane", "Methanol"}) {
        EXPECT_EQ(hybridization_codes(name), compound_hybridizations(name)) << name;
    }
}

TEST(CarbonHybridization, NonCarbonAtomsHaveNone) {
    const auto molecule = parse_ok(compound_smiles("Acetonitrile"));
    ASSERT_EQ(molecule.atoms().size(), 3U);
    EXPECT_EQ(molecule.atoms()[2].element, Element::N);
    EXPECT_EQ(carbon_hybridization(molecule, 2), std::nullopt);
}

TEST(CarbonHybridization, ChargedCarbonsAreLeftUnassigned) {
    // The bond rule disagrees with RDKit only on charged carbons (carbocations, ylide
    // carbanions), so ARISE does not assign them (they may count as any hybridization).
    // Methyl isocyanide's carbon carries a -1 charge.
    const auto molecule = parse_ok(compound_smiles("Methyl isocyanide"));
    ASSERT_EQ(molecule.atoms()[0].element, Element::C);
    ASSERT_NE(molecule.atoms()[0].charge, 0);
    EXPECT_EQ(carbon_hybridization(molecule, 0), std::nullopt);
    EXPECT_EQ(carbon_hybridization(molecule, 2), Hybridization::Sp3);  // the methyl
}

// RDKit's codes for the test compound `name` (atoms of `molecule`), with charged carbons
// set to 0, since ARISE leaves them unassigned.
std::vector<int> rdkit_codes_for_uncharged(std::string_view name, const Molecule& molecule) {
    auto codes = compound_hybridizations(name);
    for (std::size_t i = 0; i < codes.size() && i < molecule.atoms().size(); ++i) {
        if (molecule.atoms()[i].charge != 0) {
            codes[i] = 0;
        }
    }
    return codes;
}

TEST(CarbonHybridization, MatchesRdkitOnEveryUnchargedCarbonOfTheTestCompounds) {
    std::size_t carbons = 0;
    for (const auto& name : arise::test::compound_names()) {
        const auto parsed = arise::chem::parse_smiles(compound_smiles(name));
        if (!parsed.ok()) {
            continue;  // mixtures: the reader's own tests cover them
        }
        const auto& molecule = parsed.value();
        std::vector<int> ours;
        for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
            ours.push_back(hybridization_code(molecule, i));
        }
        carbons += static_cast<std::size_t>(
            std::count_if(ours.begin(), ours.end(), [](int code) { return code != 0; }));
        EXPECT_EQ(ours, rdkit_codes_for_uncharged(name, molecule)) << name;
    }
    EXPECT_GT(carbons, 500U);
}

// ---- Aromatic rings --------------------------------------------------------------

TEST(AromaticRings, FindsKekuleBenzene) {
    const auto rings =
        aromatic_rings(parse_ok(compound_smiles("Toluene (atom order 3)")));  // toluene
    ASSERT_EQ(rings.size(), 1U);
    EXPECT_EQ(rings[0], (std::vector<std::size_t>{1, 2, 3, 4, 5, 6}));
}

TEST(AromaticRings, FindsPyridine) {
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Pyridine"))).size(), 1U);
}

TEST(AromaticRings, FindsBothRingsOfKekuleNaphthalene) {
    // Written this way, the right-hand ring has only two double bonds of its own; the
    // third is the fusion bond's partner in the left ring. Both rings are aromatic.
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Naphthalene"))).size(), 2U);
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Biphenyl"))).size(), 2U);  // biphenyl
}

TEST(AromaticRings, FindsFiveMemberedHeteroaromatics) {
    // Huckel 6 pi: two ring C=C plus the lone pair of O, S or N(H/R) (user decision:
    // furan, pyrrole and thiophene are aromatic).
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Furan"))).size(), 1U);    // furan
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Pyrrole"))).size(), 1U);  // pyrrole
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("1-Methylpyrrole"))).size(),
              1U);  // N-methylpyrrole
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Thiophene"))).size(), 1U);  // thiophene
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Furan")))[0].size(), 5U);
    // Indole: the benzene ring and the pyrrole ring, which shares a fusion bond. In both
    // Kekule forms the fusion carbons' double bonds lie in the benzene ring.
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Indole"))).size(), 2U);
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Indole (atom order 2)"))).size(), 2U);
}

TEST(AromaticRings, IgnoresNonAromaticFiveRings) {
    EXPECT_TRUE(
        aromatic_rings(parse_ok(compound_smiles("2,3-Dihydrofuran"))).empty());  // 2,3-dihydrofuran
    EXPECT_TRUE(aromatic_rings(parse_ok(compound_smiles("Cyclopentadiene")))
                    .empty());  // cyclopentadiene: no lone pair
    EXPECT_TRUE(aromatic_rings(parse_ok(compound_smiles("Tetrahydrofuran"))).empty());  // THF
    // A ring S with two extra bonds has no free lone pair for the ring.
    EXPECT_TRUE(aromatic_rings(parse_ok(compound_smiles("a'-Dimethyl thiophene"))).empty());
    // Fluoranthene: three benzenoid rings; its five-membered ring has five carbons and no
    // lone-pair donor, so it is not counted.
    EXPECT_EQ(aromatic_rings(parse_ok(compound_smiles("Fluoranthene"))).size(), 3U);
}

TEST(AromaticRings, IgnoresNonAromaticSixRings) {
    EXPECT_TRUE(aromatic_rings(parse_ok(compound_smiles("Cyclohexene"))).empty());  // cyclohexene
    EXPECT_TRUE(
        aromatic_rings(parse_ok(compound_smiles("1,3-Cyclohexadiene"))).empty());  // cyclohexadiene
    EXPECT_TRUE(aromatic_rings(parse_ok(compound_smiles("Cyclohexane"))).empty());  // cyclohexane
    EXPECT_TRUE(aromatic_rings(parse_ok(compound_smiles("Quinone"))).empty());      // benzoquinone
    EXPECT_TRUE(
        aromatic_rings(parse_ok(compound_smiles("Schembl4657715"))).empty());  // thiabenzene ylide
}

TEST(AromaticRings, DoNotDependOnHowTheSmilesIsWritten) {
    const auto a = aromatic_rings(parse_ok(compound_smiles("Toluene (atom order 2)")));
    const auto b = aromatic_rings(parse_ok(compound_smiles("Toluene")));
    ASSERT_EQ(a.size(), 1U);
    ASSERT_EQ(b.size(), 1U);
    EXPECT_EQ(a[0], b[0]);  // both rings are atoms 0-5
}

// ---- Symmetry classes ------------------------------------------------------------

TEST(Symmetry, NumbersClassesInOrderOfFirstAppearance) {
    // Toluene: CH3, ipso, ortho, meta, para, meta, ortho.
    EXPECT_EQ(symmetry_classes(parse_ok(compound_smiles("Toluene (atom order 3)"))),
              (std::vector<std::size_t>{0, 1, 2, 3, 4, 3, 2}));
}

TEST(Symmetry, KekuleBenzeneStaysSymmetric) {
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Benzene")), (std::vector<int>{6}));
    // o-Xylene: the two ortho CH and the two ipso C must pair up although one ipso
    // carbon has a written double bond to its neighbour and the other a single bond.
    EXPECT_EQ(carbon_class_sizes(compound_smiles("O-Xylene")), (std::vector<int>{2, 2, 2, 2}));
}

TEST(Symmetry, RegioisomersGiveTheirTextbookCarbonCounts) {
    // Distinct 13C environments: toluene 5, o-/m-/p-xylene 4/5/3, ethylbenzene 6.
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Toluene (atom order 3)")).size(), 5U);
    EXPECT_EQ(carbon_class_sizes(compound_smiles("O-Xylene")).size(), 4U);
    EXPECT_EQ(carbon_class_sizes(compound_smiles("M-Xylene")), (std::vector<int>{1, 1, 2, 2, 2}));
    EXPECT_EQ(carbon_class_sizes(compound_smiles("P-Xylene (atom order 2)")),
              (std::vector<int>{2, 2, 4}));
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Ethylbenzene")).size(), 6U);
}

TEST(Symmetry, AliphaticEquivalences) {
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Isopropanol")),
              (std::vector<int>{1, 2}));  // isopropanol
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Diethyl Ether")),
              (std::vector<int>{2, 2}));  // diethyl ether
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Neopentane")),
              (std::vector<int>{1, 4}));  // neopentane
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Ethyl Acetate")), (std::vector<int>{1, 1, 1, 1}));
}

TEST(Symmetry, FiveMemberedHeteroaromaticsStaySymmetric) {
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Furan")), (std::vector<int>{2, 2}));  // furan
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Thiophene")),
              (std::vector<int>{2, 2}));  // thiophene
    EXPECT_EQ(carbon_class_sizes(compound_smiles("2,5-Dimethylfuran")),
              (std::vector<int>{2, 2, 2}));  // 2,5-dimethylfuran
}

TEST(Symmetry, NaphthaleneHasThreeCarbonEnvironments) {
    // alpha (4), beta (4) and the two fusion carbons.
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Naphthalene")), (std::vector<int>{2, 4, 4}));
}

TEST(Symmetry, DoesNotDependOnHowTheSmilesIsWritten) {
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Toluene (atom order 3)")),
              carbon_class_sizes(compound_smiles("Toluene (atom order 2)")));
    EXPECT_EQ(carbon_class_sizes(compound_smiles("P-Xylene (atom order 2)")),
              carbon_class_sizes(compound_smiles("P-Xylene")));
}

TEST(Symmetry, ExplicitHydrogensDoNotChangeEquivalence) {
    // Ethane with one methyl's hydrogens written as atoms: both CH3 are still equivalent.
    EXPECT_EQ(carbon_class_sizes(syntax_case("explicit hydrogen atoms on one ethane methyl")),
              (std::vector<int>{2}));
}

TEST(Symmetry, HydrogenCountDistinguishesOtherwiseAlikeAtoms) {
    // Propane with a CH2 radical end: the two end carbons differ only in hydrogens.
    EXPECT_EQ(carbon_class_sizes(syntax_case("CH2 radical end of propane")),
              (std::vector<int>{1, 1, 1}));
}

TEST(Symmetry, ResonanceSharedOxygensAreEquivalent) {
    // The two O of a nitro group or a carboxylate share the charge by resonance.
    const auto nitro = parse_ok(compound_smiles("Nitroethane"));
    const auto classes = arise::chem::symmetry_classes(nitro);
    EXPECT_EQ(classes[3], classes[4]);
    const auto glycine = parse_ok(compound_smiles("alpha-Glycine"));
    const auto zwitterion = arise::chem::symmetry_classes(glycine);
    EXPECT_EQ(zwitterion[2], zwitterion[3]);
    // ...but an acid's OH and C=O stay distinct.
    const auto acid = arise::chem::symmetry_classes(parse_ok(compound_smiles("Acetic Acid")));
    EXPECT_NE(acid[2], acid[3]);
}

TEST(Symmetry, IsotopeLabelBreaksEquivalence) {
    // 13C-labelled methyl in propane: the two methyls are no longer equivalent.
    EXPECT_EQ(carbon_class_sizes(compound_smiles("Propane (1-13C)")), (std::vector<int>{1, 1, 1}));
}

// ---- NMR equivalence ---------------------------------------------------------------
// Expected 13C environment counts under the chemist's rules (perception.hpp,
// NmrEquivalence). The rules themselves are checked against nmrshiftdb2's assigned spectra
// in reference_test.cpp.

struct Range {
    int fewest;
    int most;
    friend bool operator==(const Range&, const Range&) = default;
};

std::ostream& operator<<(std::ostream& out, const Range& range) {
    return out << "[" << range.fewest << ", " << range.most << "]";
}

Range environments(std::string_view smiles) {
    const auto molecule = parse_ok(smiles);
    const auto equivalence = arise::chem::nmr_equivalence(molecule);
    return Range{
        .fewest = arise::chem::carbon_environments(molecule, equivalence.fewest),
        .most = arise::chem::carbon_environments(molecule, equivalence.most),
    };
}

int topological(std::string_view smiles) {
    const auto molecule = parse_ok(smiles);
    return arise::chem::carbon_environments(molecule, symmetry_classes(molecule));
}

TEST(NmrEquivalence, EZGeometrySplitsAlways) {
    EXPECT_EQ(topological(compound_smiles("Linalool")), 9);  // linalool
    EXPECT_EQ(environments(compound_smiles("Linalool")), (Range{10, 10}));
    EXPECT_EQ(environments(compound_smiles("2-Methyl-2-butene")),
              (Range{5, 5}));  // 2-methyl-2-butene
    EXPECT_EQ(environments(compound_smiles("Cyclohexanone Oxime")),
              (Range{6, 6}));  // cyclohexanone oxime
    EXPECT_EQ(environments(compound_smiles("Acetone, oxime")), (Range{3, 3}));  // acetone oxime
}

TEST(NmrEquivalence, EZNeedsDifferentGroupsOnTheOtherEnd) {
    EXPECT_EQ(environments(compound_smiles("2,3-Dimethyl-2-butene")),
              (Range{2, 2}));  // tetramethylethylene
    EXPECT_EQ(environments(compound_smiles("Isobutylene")),
              (Range{3, 3}));  // isobutylene: =CH2 end
}

TEST(NmrEquivalence, AmideNGroupsAlwaysSplit) {
    EXPECT_EQ(topological(compound_smiles("Dimethylformamide")), 2);  // DMF
    EXPECT_EQ(environments(compound_smiles("Dimethylformamide")), (Range{3, 3}));
    EXPECT_EQ(environments(compound_smiles("Dimethylacetamide")),
              (Range{4, 4}));  // N,N-dimethylacetamide
    EXPECT_EQ(environments(compound_smiles("N,N-Dimethylthioformamide")),
              (Range{3, 3}));  // thioamide
}

TEST(NmrEquivalence, CarbamateUreaAndEnaminoneNGroupsGiveARange) {
    EXPECT_EQ(environments(compound_smiles("Tetramethylurea")), (Range{2, 3}));  // tetramethylurea
    EXPECT_EQ(environments(compound_smiles("Carbamic acid, dimethyl-, methyl ester")),
              (Range{3, 4}));  // methyl dimethylcarbamate
    EXPECT_EQ(environments(compound_smiles("4-(Dimethylamino)but-3-en-2-one")),
              (Range{5, 6}));  // enaminone
}

// The chemist's decision of 2026-11-08: splits that nmrshiftdb2's assigned spectra
// support only weakly (Reference.AlwaysSplitCarbonsAreSeparatePeaksWhereNmrshiftdb2ListsBoth
// counts them) give a range instead of always splitting.
TEST(NmrEquivalence, WeaklySupportedSplitsGiveARange) {
    // Amidine/guanidine C=N: the two N(CH3)2 groups may or may not differ.
    EXPECT_EQ(environments(compound_smiles("1,1,3,3-Tetramethylguanidine")), (Range{2, 3}));
    // C=C whose far end carries an amino N (push-pull): the two acetyls may differ, and
    // the N-methyls (an enaminone N) may differ too.
    EXPECT_EQ(environments(compound_smiles("3-[(Dimethylamino)methylidene]pentane-2,4-dione")),
              (Range{5, 8}));
    // Push-pull C=C with the amino groups on the identical end (a ketene aminal): the two
    // N(CH3)2 groups may differ, and so may the methyls within each.
    EXPECT_EQ(environments(compound_smiles("3-Buten-2-one, 4,4-bis(dimethylamino)-")),
              (Range{5, 8}));
    // N-acyl azole: the two ring halves of N-acetylpyrrole may differ.
    EXPECT_EQ(environments(compound_smiles("N-Acetylpyrrole")), (Range{4, 6}));
    // N-acyl azole that ARISE does not perceive as aromatic (a 1,2,3-triazole).
    EXPECT_EQ(environments(compound_smiles("1-(2H-1,2,3-Triazol-2-yl)ethanone")), (Range{3, 4}));
}

TEST(NmrEquivalence, DiazoCarbonNeverSplitsItsGroups) {
    // C=N+=N-: nmrshiftdb2 lists the two tert-butyl groups at one shift in every case.
    EXPECT_EQ(environments(compound_smiles("3-Diazo-2,2,4,4-tetramethylpentane")), (Range{3, 3}));
}

// The chemist's decisions of 2026-11-09, from nmrshiftdb2's assigned spectra (each pattern's
// support is in the plan; the always-split one is checked by
// Reference.AlwaysSplitCarbonsAreSeparatePeaksWhereNmrshiftdb2ListsBoth).
TEST(NmrEquivalence, NitrosamineNGroupsAlwaysSplit) {
    EXPECT_EQ(environments(compound_smiles("N-Nitrosodimethylamine")),
              (Range{2, 2}));  // N-nitrosodimethylamine
    EXPECT_EQ(environments(compound_smiles("N-Nitrosopiperidine")),
              (Range{5, 5}));  // N-nitrosopiperidine
}

// The chemist's decisions of 2026-11-11: three existing rules extended to forms they missed
// (records 20048921, 90290 and 20025993 of nmrshiftdb2 list these splits).
TEST(NmrEquivalence, ChargeSeparatedNitrosamineNGroupsAlwaysSplit) {
    EXPECT_EQ(environments(compound_smiles("N-Nitrosodimethylamine (charge-separated)")),
              (Range{2, 2}));
}

TEST(NmrEquivalence, AnAxiallyChiralAlleneCountsAsAStereocentre) {
    // Both allene terminals carry two different groups, so the isopropyl methyls may differ.
    EXPECT_EQ(environments(compound_smiles("Methyl 5-methylhexa-2,3-dienoate")), (Range{7, 8}));
    // One terminal with two identical groups: no chirality element from the allene, so the
    // isopropyl methyls stay together (only the terminal pair may split: 6 to 7, not 8).
    EXPECT_EQ(environments(compound_smiles("Schembl4340083")), (Range{6, 7}));
}

TEST(NmrEquivalence, ASulfilimineSulfurIsAStereocentre) {
    // S(=N) with two different groups, like a sulfoxide S(=O): the isopropyl methyls may differ.
    EXPECT_EQ(environments(compound_smiles(
                  "tert-butyl N-[methyl(propan-2-yl)-lambda4-sulfanylidene]carbamate")),
              (Range{6, 7}));
}

TEST(NmrEquivalence, GemPairOnARingAtomGivesARange) {
    // 1,1-dimethylcyclohexane: the two exocyclic methyls may differ.
    EXPECT_EQ(environments(compound_smiles("1,1-Dimethylcyclohexane")), (Range{5, 6}));
    // A ring atom's two ring neighbours are not a gem pair: cyclohexanol keeps its halves.
    EXPECT_EQ(environments(compound_smiles("Cyclohexanol")), (Range{4, 4}));
    // An open-chain atom is not a ring atom: propane's two methyls stay one environment.
    EXPECT_EQ(environments(compound_smiles("Propane")), (Range{2, 2}));
}

TEST(NmrEquivalence, AlleneTerminalPairGivesARange) {
    // 2-methylpenta-2,3-diene: the two methyls on one terminal may differ.
    EXPECT_EQ(environments(compound_smiles("2,3-Pentadiene, 2-methyl-")), (Range{5, 6}));
}

TEST(NmrEquivalence, DiastereotopicGroupsGiveARange) {
    EXPECT_EQ(environments(compound_smiles("DL-valine")), (Range{4, 5}));  // valine
    EXPECT_EQ(environments(compound_smiles("4-Methyl-2-Pentanol")),
              (Range{5, 6}));  // 4-methyl-2-pentanol
    EXPECT_EQ(environments(compound_smiles("Propane, 2-(ethenylsulfinyl)-")),
              (Range{4, 5}));  // chiral sulfoxide
    EXPECT_EQ(environments(compound_smiles("[Methyl(propan-2-yl)phosphoryl]benzene")),
              (Range{7, 8}));  // chiral P
}

TEST(NmrEquivalence, SeveralProchiralCentresAlsoGiveARange) {
    // nmrshiftdb2 10017805: no stereocentre, yet each acetal's two OEt groups give
    // different peaks, because several prochiral centres make each other's groups
    // diastereotopic. Over-splitting only widens the range, which is the safe direction.
    const auto range =
        environments(compound_smiles("2-(Diethoxymethyl)-1,1-diethoxy-3-methylbutane"));
    EXPECT_EQ(range.fewest, 6);
    EXPECT_GE(range.most, 8);
}

TEST(NmrEquivalence, OneProchiralCentreAloneKeepsEquivalence) {
    // Isopropanol: one prochiral centre and no stereocentre, so no rule splits its methyls.
    EXPECT_EQ(environments(compound_smiles("Isopropanol")), (Range{2, 2}));
}

TEST(NmrEquivalence, UnknownStereoCanBreakASymmetry) {
    // A hexitol written without stereo (mannitol and sorbitol share this SMILES): every
    // stereocentre is split on its own in the range, so the range spans the topological
    // count and full splitting.
    EXPECT_EQ(environments(compound_smiles("Hexitol")), (Range{3, 6}));
}

TEST(NmrEquivalence, NoRuleSplitsTertButylMethylsOrPhenylOrthoCarbons) {
    // tert-butyl methyls (three identical groups) and a phenyl's ortho carbons: no rule
    // splits them, even next to a stereocentre.
    EXPECT_EQ(environments(compound_smiles("3,3-Dimethyl-2-butanol")),
              (Range{4, 4}));  // 3,3-dimethyl-2-butanol
    EXPECT_EQ(environments(compound_smiles("1-Phenylethanol")), (Range{6, 6}));  // 1-phenylethanol
    EXPECT_EQ(environments(compound_smiles("Ethyl Acetate")), (Range{4, 4}));    // ethyl acetate
    EXPECT_EQ(environments(compound_smiles("P-Xylene (atom order 2)")), (Range{3, 3}));  // p-xylene
}

TEST(NmrEquivalence, ProchiralCentresBeyondCarbonMaySplit) {
    // A quaternary N+ with two methyls and a stereocentre elsewhere: the "most" partition
    // splits the methyls. Lauryl hydroxysultaine (PubChem CID 114474; the substitute,
    // 2026-11-21, for a hand-made zwitterion with no PubChem entry) has 17 carbons: 12 in
    // the lauryl chain, the two N-methyls (one topological class) and 3 in the
    // hydroxypropyl sulfonate, so 16..17 carbon environments.
    EXPECT_EQ(environments(compound_smiles("Lauryl Hydroxysultaine")), (Range{16, 17}));
    // A neutral amine inverts fast, so its two methyls stay equivalent.
    EXPECT_EQ(environments(compound_smiles("Dimepranol")), (Range{4, 4}));
}

TEST(NmrEquivalence, NeverFewerClassesThanTopology) {
    for (const auto& smiles : {
             compound_smiles("DL-valine"),
             compound_smiles("Dimethylformamide"),
             compound_smiles("Hexitol"),
         }) {
        const auto range = environments(smiles);
        EXPECT_GE(range.fewest, topological(smiles)) << smiles;
        EXPECT_GE(range.most, range.fewest) << smiles;
    }
}

// ---- Functional groups (docs/chem_group_definitions.md) -----------------------------
// One test per family; each checks the inclusions and exclusions the definitions table
// names. External checks against ChEBI and RDKit live in reference_test.cpp.

TEST(Groups, AcidsAldehydesKetones) {
    EXPECT_EQ(groups(compound_smiles("Acetic Acid"), Group::CarboxylicAcid), 1);
    EXPECT_EQ(groups(compound_smiles("Formic Acid"), Group::CarboxylicAcid), 1);  // formic acid...
    EXPECT_EQ(groups(compound_smiles("Formic Acid"), Group::Aldehyde), 0);  // ...not an aldehyde
    EXPECT_EQ(groups(compound_smiles("Benzoic Acid"), Group::CarboxylicAcid), 1);  // benzoic acid
    EXPECT_EQ(groups(compound_smiles("Acetaldehyde"), Group::Aldehyde), 1);
    EXPECT_EQ(groups(compound_smiles("Formaldehyde"), Group::Aldehyde),
              1);  // formaldehyde (user decision)
    EXPECT_EQ(groups(compound_smiles("Methyl formate"), Group::Aldehyde), 0);     // methyl formate
    EXPECT_EQ(groups(compound_smiles("Dimethylformamide"), Group::Aldehyde), 0);  // DMF
    EXPECT_EQ(groups(compound_smiles("Acetone"), Group::Ketone), 1);
    EXPECT_EQ(groups(compound_smiles("Acetophenone"), Group::Ketone), 1);  // acetophenone
    EXPECT_EQ(groups(compound_smiles("Quinone"), Group::Ketone), 2);       // benzoquinone
    EXPECT_EQ(groups(compound_smiles("Ethyl Acetate"), Group::Ketone), 0);
    EXPECT_EQ(groups(compound_smiles("Acetaldehyde"), Group::Ketone),
              0);  // an aldehyde has one carbon side
}

TEST(Groups, EstersAndTheirNeighbours) {
    EXPECT_EQ(groups(compound_smiles("Ethyl Acetate"), Group::CarboxylicEster),
              1);  // ethyl acetate
    EXPECT_EQ(groups(compound_smiles("Methyl formate"), Group::CarboxylicEster),
              1);  // formate (user decision)
    EXPECT_EQ(groups(compound_smiles("Gamma-Butyrolactone"), Group::CarboxylicEster),
              1);  // lactone
    EXPECT_EQ(groups(compound_smiles("Dimethyl Carbonate"), Group::CarboxylicEster),
              0);                                                                       // carbonate
    EXPECT_EQ(groups(compound_smiles("Methyl Carbamate"), Group::CarboxylicEster), 0);  // carbamate
    EXPECT_EQ(groups(compound_smiles("Acetic Anhydride"), Group::CarboxylicEster), 0);  // anhydride
    EXPECT_EQ(groups(compound_smiles("Dimethyl Carbonate"), Group::CarbonateEster), 1);
    EXPECT_EQ(groups(compound_smiles("Ethylene carbonate"), Group::CarbonateEster),
              1);  // ethylene carbonate
    EXPECT_EQ(groups(compound_smiles("Ethyl Acetate"), Group::CarbonateEster), 0);
    EXPECT_EQ(groups(compound_smiles("Methyl hydrogen carbonate"), Group::CarbonateEster),
              0);  // monoester: one O is an OH
    EXPECT_EQ(groups(compound_smiles("Acetic Anhydride"), Group::CarboxylicAnhydride), 1);
    EXPECT_EQ(groups(compound_smiles("Succinic anhydride"), Group::CarboxylicAnhydride),
              1);  // succinic anhydride
    EXPECT_EQ(groups(compound_smiles("Ethyl Acetate"), Group::CarboxylicAnhydride), 0);
}

TEST(Groups, AmidesCarbamatesUreas) {
    EXPECT_EQ(groups(compound_smiles("Acetamide"), Group::Carboxamide), 1);
    EXPECT_EQ(groups(compound_smiles("Dimethylformamide"), Group::Carboxamide), 1);  // formamide
    EXPECT_EQ(groups(compound_smiles("N-Methylacetamide"), Group::Carboxamide), 1);
    EXPECT_EQ(groups(compound_smiles("2-Pyrrolidone"), Group::Carboxamide), 1);  // lactam
    EXPECT_EQ(groups(compound_smiles("Urea"), Group::Carboxamide), 0);  // urea (user decision)
    EXPECT_EQ(groups(compound_smiles("Methyl Carbamate"), Group::Carboxamide),
              0);  // carbamate (user decision)
    EXPECT_EQ(groups(compound_smiles("Methyl Carbamate"), Group::CarbamateEster), 1);
    EXPECT_EQ(groups(compound_smiles("2-Oxazolidinone"), Group::CarbamateEster),
              1);  // oxazolidinone
    EXPECT_EQ(groups(compound_smiles("Urea"), Group::CarbamateEster), 0);
    EXPECT_EQ(groups(compound_smiles("Carbamic Acid"), Group::CarbamateEster),
              0);  // carbamic acid: no ester O
    // O-carbamoyl oximes count (user decision 2026-10-31, as in ChEBI): methomyl, aldicarb.
    EXPECT_EQ(groups(compound_smiles("Methomyl"), Group::CarbamateEster), 1);
    EXPECT_EQ(groups(compound_smiles("Aldicarb"), Group::CarbamateEster), 1);
    EXPECT_EQ(groups(compound_smiles("Methomyl"), Group::CarboxylicEster), 0);
    EXPECT_EQ(groups(compound_smiles("Methylamino methylcarbamate"), Group::CarbamateEster),
              0);  // O on an amine N: not oxime
    EXPECT_EQ(groups(compound_smiles("propan-2-one O-acetyl oxime"), Group::CarboxylicEster),
              0);  // O-acyl oxime: not decided
    EXPECT_EQ(groups(compound_smiles("Urea"), Group::Urea), 1);
    EXPECT_EQ(groups(compound_smiles("Tetramethylurea"), Group::Urea), 1);
    EXPECT_EQ(groups(compound_smiles("2-Imidazolidinone"), Group::Urea), 1);  // imidazolidinone
    EXPECT_EQ(groups(compound_smiles("Acetamide"), Group::Urea), 0);
}

TEST(Groups, HydroxylFamilies) {
    EXPECT_EQ(groups(compound_smiles("Ethanol"), Group::Alcohol), 1);
    EXPECT_EQ(groups(compound_smiles("Vinyl alcohol"), Group::Alcohol), 0);
    EXPECT_EQ(groups(compound_smiles("Phenol"), Group::Alcohol), 0);
    EXPECT_EQ(groups(compound_smiles("Acetic Acid"), Group::Alcohol), 0);
    EXPECT_EQ(groups(compound_smiles("Vinyl alcohol"), Group::Enol), 1);
    EXPECT_EQ(groups(compound_smiles("Ethanol"), Group::Enol), 0);
    EXPECT_EQ(groups(compound_smiles("Phenol"), Group::Enol), 0);
    EXPECT_EQ(groups(compound_smiles("Methylketene"), Group::Enol),
              0);  // ketene: the O carries no H
    EXPECT_EQ(groups(compound_smiles("Phenol"), Group::Phenol), 1);
    EXPECT_EQ(groups(compound_smiles("2-Naphthol"), Group::Phenol), 1);  // 2-naphthol
    EXPECT_EQ(groups(compound_smiles("4-Hydroxybenzoic Acid"), Group::Phenol), 1);
    EXPECT_EQ(groups(compound_smiles("3-Hydroxypyridine"), Group::Phenol), 0);  // hydroxypyridine
    EXPECT_EQ(groups(compound_smiles("3-Hydroxypyridine"), Group::Pyridone), 1);
    EXPECT_EQ(groups(compound_smiles("2-Pyridone"), Group::Pyridone), 1);         // 2-pyridone
    EXPECT_EQ(groups(compound_smiles("4-Hydroxypyridine"), Group::Pyridone), 1);  // 4-pyridone
    EXPECT_EQ(groups(compound_smiles("Phenol"), Group::Pyridone), 0);
    EXPECT_EQ(groups(compound_smiles("2-Piperidone"), Group::Pyridone),
              0);  // valerolactam: no ring C=C
}

TEST(Groups, Ethers) {
    EXPECT_EQ(groups(compound_smiles("Diethyl Ether"), Group::Ether), 1);
    EXPECT_EQ(groups(compound_smiles("Methoxyacetone"), Group::Ether), 1);   // methoxyacetone
    EXPECT_EQ(groups(compound_smiles("Anisole"), Group::Ether), 1);          // anisole
    EXPECT_EQ(groups(compound_smiles("1,4-Dioxane"), Group::Ether), 2);      // dioxane
    EXPECT_EQ(groups(compound_smiles("Tetrahydrofuran"), Group::Ether), 1);  // THF
    EXPECT_EQ(groups(compound_smiles("Propylene Oxide"), Group::Ether), 1);  // epoxide
    EXPECT_EQ(groups(compound_smiles("Ethyl Acetate"), Group::Ether), 0);    // ester O
    EXPECT_EQ(groups(compound_smiles("Furan"), Group::Ether), 0);  // furan O (user decision)
}

// Methoxy (the chemist's decisions, 2026-11-22): a CH3 single-bonded to an O that is bonded
// to exactly one other carbon (aliphatic, aryl or carbonyl), counted per OCH3; it matches
// RDKit's fr_methoxy, which also counts esters and aryl OMe.
TEST(Groups, Methoxy) {
    EXPECT_EQ(groups(compound_smiles("Anisole"), Group::Methoxy), 1);         // aryl OMe
    EXPECT_EQ(groups(compound_smiles("Methoxyacetone"), Group::Methoxy), 1);  // aliphatic OMe
    EXPECT_EQ(groups(compound_smiles("Methyl formate"), Group::Methoxy), 1);  // methyl ester
    EXPECT_EQ(groups(compound_smiles("Dimethyl Carbonate"), Group::Methoxy), 2);
    EXPECT_EQ(groups(compound_smiles("Methyl Carbamate"), Group::Methoxy), 1);
    EXPECT_EQ(groups(compound_smiles("Dimethyl Ether"), Group::Methoxy), 2);  // per OCH3
    EXPECT_EQ(groups(compound_smiles("Methanol"), Group::Methoxy), 0);  // O has no other carbon
    EXPECT_EQ(groups(compound_smiles("Diethyl Ether"), Group::Methoxy), 0);  // no CH3 on O
    EXPECT_EQ(groups(compound_smiles("Ethyl Acetate"), Group::Methoxy), 0);  // ethyl ester
    EXPECT_EQ(groups(compound_smiles("N,O-Dimethylhydroxylamine"), Group::Methoxy),
              0);  // the O's other atom is N
    EXPECT_EQ(groups(compound_smiles("Silane, methoxytrimethyl-"), Group::Methoxy),
              0);  // the O's other atom is Si
    EXPECT_EQ(groups(compound_smiles("Thioanisole"), Group::Methoxy), 0);  // S, not O
    EXPECT_EQ(groups(compound_smiles("Silane, methoxy-"), Group::Methoxy),
              0);  // CH3-O-SiH3: the O's other atom is Si, and the SiH3 is not a methyl
}

TEST(Groups, NitrogenFamilies) {
    EXPECT_EQ(groups(compound_smiles("Methylamine"), Group::Amine), 1);
    EXPECT_EQ(groups(compound_smiles("Triethylamine"), Group::Amine), 1);
    EXPECT_EQ(groups(compound_smiles("Aniline"), Group::Amine), 1);  // aniline
    EXPECT_EQ(groups(compound_smiles("Acetamide"), Group::Amine), 0);
    EXPECT_EQ(groups(compound_smiles("Pyridine"), Group::Amine), 0);  // pyridine N
    EXPECT_EQ(groups(compound_smiles("Pyrrole"), Group::Amine), 0);   // pyrrole N (user decision)
    EXPECT_EQ(groups(compound_smiles("Methanesulfonamide"), Group::Amine), 0);   // sulfonamide N
    EXPECT_EQ(groups(compound_smiles("Methylhydrazine"), Group::Amine), 0);      // hydrazine N
    EXPECT_EQ(groups(compound_smiles("Methylhydroxylamine"), Group::Amine), 0);  // hydroxylamine N
    EXPECT_EQ(groups(compound_smiles("Thioacetamide"), Group::Amine), 0);        // thioamide N
    EXPECT_EQ(groups(compound_smiles("Methanesulfonamide"), Group::Sulfonamide), 1);
    EXPECT_EQ(groups(compound_smiles("N,N-Dimethylmethanesulphonamide"), Group::Sulfonamide), 1);
    EXPECT_EQ(groups(compound_smiles("Methanesulfinamide"), Group::Sulfonamide),
              0);  // sulfinamide: one S=O
    EXPECT_EQ(groups(compound_smiles("Methylhydrazine"), Group::Hydrazine), 1);
    EXPECT_EQ(groups(compound_smiles("Acetohydrazide"), Group::Hydrazine), 1);  // hydrazide
    EXPECT_EQ(groups(compound_smiles("Methylhydroxylamine"), Group::Hydroxylamine), 1);
    EXPECT_EQ(groups(compound_smiles("Acetone, oxime"), Group::Hydroxylamine), 0);  // oxime
    EXPECT_EQ(groups(compound_smiles("2-Propanamine, N-methyl-"), Group::Imine), 1);
    EXPECT_EQ(groups(compound_smiles("Acetone, oxime"), Group::Imine), 0);
    EXPECT_EQ(groups(compound_smiles("2-Propanone, hydrazone"), Group::Imine), 0);  // hydrazone
    EXPECT_EQ(groups(compound_smiles("Acetone, oxime"), Group::Oxime), 1);
    EXPECT_EQ(groups(compound_smiles("2-Propanamine, N-methyl-"), Group::Oxime), 0);
    EXPECT_EQ(groups(compound_smiles("2-Propanone, o-methyloxime"), Group::Oxime),
              0);  // O-methyl oxime: no OH
    EXPECT_EQ(groups(compound_smiles("Acetonitrile"), Group::Nitrile), 1);
}

TEST(Groups, AlkenesAndAromaticRings) {
    EXPECT_EQ(groups(compound_smiles("1,3-Butadiene"), Group::Alkene), 2);
    EXPECT_EQ(groups(compound_smiles("Styrene"), Group::Alkene),
              1);  // styrene: ring bonds excluded
    EXPECT_EQ(groups(compound_smiles("Benzene"), Group::Alkene), 0);
    EXPECT_EQ(groups(compound_smiles("Quinone"), Group::Alkene), 2);  // quinone is not aromatic
    EXPECT_EQ(groups(compound_smiles("Furan"), Group::Alkene), 0);    // furan ring bonds
    EXPECT_EQ(groups(compound_smiles("2-Ethenylfuran"), Group::Alkene), 1);  // 2-vinylfuran
    EXPECT_EQ(groups(compound_smiles("Styrene"), Group::BenzeneRing), 1);
    EXPECT_EQ(groups(compound_smiles("Naphthalene"), Group::BenzeneRing), 2);
    EXPECT_EQ(groups(compound_smiles("Pyridine"), Group::BenzeneRing), 0);
    EXPECT_EQ(groups(compound_smiles("Pyridine"), Group::PyridineRing), 1);
    EXPECT_EQ(groups(compound_smiles("Quinoline"), Group::PyridineRing), 1);  // quinoline
    EXPECT_EQ(groups(compound_smiles("Quinoline"), Group::BenzeneRing), 1);
    EXPECT_EQ(groups(compound_smiles("Furan"), Group::FuranRing), 1);
    EXPECT_EQ(groups(compound_smiles("Pyrrole"), Group::FuranRing),
              0);  // the heteroatom decides the ring
    EXPECT_EQ(groups(compound_smiles("Thiophene"), Group::PyrroleRing), 0);
    EXPECT_EQ(groups(compound_smiles("Benzofuran"), Group::FuranRing), 1);       // benzofuran
    EXPECT_EQ(groups(compound_smiles("Tetrahydrofuran"), Group::FuranRing), 0);  // THF
    EXPECT_EQ(groups(compound_smiles("Pyrrole"), Group::PyrroleRing), 1);
    EXPECT_EQ(groups(compound_smiles("1-Methylpyrrole"), Group::PyrroleRing), 1);
    EXPECT_EQ(groups(compound_smiles("Indole"), Group::PyrroleRing), 1);  // indole
    EXPECT_EQ(groups(compound_smiles("Indole"), Group::BenzeneRing), 1);
    EXPECT_EQ(groups(compound_smiles("Pyrrolidine"), Group::PyrroleRing), 0);  // pyrrolidine
    EXPECT_EQ(groups(compound_smiles("Thiophene"), Group::ThiopheneRing), 1);
    EXPECT_EQ(groups(compound_smiles("Tetrahydrothiophene"), Group::ThiopheneRing), 0);
}

// ---- The real PubChem C4H8O2 result ------------------------------------------------

struct Expected {
    std::string_view cid;
    std::vector<int> carbon_classes;
    std::map<Group, int> groups;  // groups not listed must be absent
};

void expect_perceived(const arise::io::Candidate& candidate, const Expected& expected) {
    const auto molecule = parse_ok(candidate.smiles);
    EXPECT_EQ(carbon_class_sizes(candidate.smiles), expected.carbon_classes) << expected.cid;
    for (const auto group : arise::chem::kAllGroups) {
        const auto it = expected.groups.find(group);
        const int want = it == expected.groups.end() ? 0 : it->second;
        EXPECT_EQ(count_group(molecule, group), want)
            << expected.cid << " " << candidate.smiles << ": " << arise::chem::to_string(group);
    }
}

TEST(Perception, ReadsKeyC4H8O2StructuresAsPubChemWritesThem) {
    // Hand-derived: carbon environments and functional groups of six C4H8O2 isomers.
    const std::vector<Expected> cases{
        {
            .cid = "CID 8857",
            .carbon_classes = {1, 1, 1, 1},
            .groups = {{Group::CarboxylicEster, 1}},
        },
        {
            .cid = "CID 11124",
            .carbon_classes = {1, 1, 1, 1},
            .groups = {{Group::CarboxylicEster, 1}, {Group::Methoxy, 1}},  // methyl ester
        },
        {
            .cid = "CID 22172",
            .carbon_classes = {1, 1, 1, 1},
            .groups = {{Group::Ketone, 1}, {Group::Ether, 1}, {Group::Methoxy, 1}},
        },
        {.cid = "CID 264", .carbon_classes = {1, 1, 1, 1}, .groups = {{Group::CarboxylicAcid, 1}}},
        {.cid = "CID 6590", .carbon_classes = {1, 1, 2}, .groups = {{Group::CarboxylicAcid, 1}}},
        {.cid = "CID 31275", .carbon_classes = {4}, .groups = {{Group::Ether, 2}}},  // dioxane
    };
    const auto imported = arise::io::import_pubchem_csv(arise::test::kPubChemC4H8O2);
    ASSERT_TRUE(imported.ok()) << imported.error().message;
    for (const auto& expected : cases) {
        const auto& candidates = imported.value().candidates;
        const auto found = std::find_if(candidates.begin(), candidates.end(),
                                        [&](const auto& c) { return c.name == expected.cid; });
        ASSERT_NE(found, candidates.end()) << expected.cid;
        expect_perceived(*found, expected);
    }
}

TEST(Perception, EveryValidC4H8O2CandidateHasFourCarbonsInItsClasses) {
    const auto imported = arise::io::import_pubchem_csv(arise::test::kPubChemC4H8O2);
    ASSERT_TRUE(imported.ok()) << imported.error().message;
    for (const auto& candidate : imported.value().candidates) {
        const auto sizes = carbon_class_sizes(candidate.smiles);
        int carbons = 0;
        for (const int size : sizes) {
            carbons += size;
        }
        EXPECT_EQ(carbons, 4) << candidate.smiles;
    }
}

TEST(Groups, ChargeSeparatedFormsCountByTheirNeutralParent) {
    // User decision 2026-11-01: carboxylate as acid, ammonium as amine, alkoxide as
    // alcohol; an N-oxide is neither amine nor hydroxylamine; an isocyanide is no nitrile;
    // a carbonyl oxide is no ketone.
    EXPECT_EQ(groups(compound_smiles("alpha-Glycine"), Group::CarboxylicAcid),
              1);  // glycine zwitterion
    EXPECT_EQ(groups(compound_smiles("alpha-Glycine"), Group::Amine), 1);
    EXPECT_EQ(groups(compound_smiles("3-Azaniumylpropan-1-olate"), Group::Alcohol), 1);
    EXPECT_EQ(groups(compound_smiles("Trimethylamine oxide"), Group::Amine),
              0);  // trimethylamine N-oxide
    EXPECT_EQ(groups(compound_smiles("Trimethylamine oxide"), Group::Hydroxylamine), 0);
    EXPECT_EQ(groups(compound_smiles("Methyl isocyanide"), Group::Nitrile),
              0);  // methyl isocyanide
    EXPECT_EQ(groups(compound_smiles("Acetonitrile"), Group::Nitrile), 1);
    EXPECT_EQ(groups(compound_smiles("Schembl6298161"), Group::Ketone), 0);  // carbonyl oxide
    EXPECT_EQ(groups(compound_smiles("Pyridine N-Oxide"), Group::PyridineRing),
              1);  // pyridine N-oxide
}

// The key atoms of each occurrence, per hallmark (checked against RDKit's fragment
// patterns by Reference.KeyAtomsMatchRdkitFragmentPatterns).
std::vector<std::vector<std::vector<std::size_t>>> key_atoms(std::string_view smiles, Group group) {
    std::vector<std::vector<std::vector<std::size_t>>> result;
    for (const auto& occurrence : arise::chem::find_group(parse_ok(smiles), group)) {
        result.push_back(occurrence.key_atoms);
    }
    return result;
}

using Keys = std::vector<std::vector<std::vector<std::size_t>>>;

TEST(Groups, KeyAtomsCarryEachHallmark) {
    // Ethyl acetate: C=O carbon 3, O-C carbon 1.
    EXPECT_EQ(key_atoms(compound_smiles("Ethyl Acetate (atom order 2)"), Group::CarboxylicEster),
              (Keys{{{3}, {1}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Acetic Anhydride"), Group::CarboxylicAnhydride),
              (Keys{{{1, 4}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Ethanol"), Group::Alcohol), (Keys{{{1}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Dimethyl Ether"), Group::Ether), (Keys{{{0, 2}}}));
    // Methoxy: one occurrence per OCH3, its CH3 carbon the key atom.
    EXPECT_EQ(key_atoms(compound_smiles("Dimethyl Ether"), Group::Methoxy), (Keys{{{0}}, {{2}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Dimethylethylamine"), Group::Amine), (Keys{{{1, 3, 4}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Acetonitrile"), Group::Nitrile), (Keys{{{1}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Propylene"), Group::Alkene), (Keys{{{1, 2}}}));
    EXPECT_EQ(key_atoms(compound_smiles("Pyridine"), Group::PyridineRing),
              (Keys{{{0, 1, 2, 4, 5}}}));
    EXPECT_EQ(key_atoms(compound_smiles("2-Pyridone"), Group::Pyridone), (Keys{{{1}}}));
    // A sulfonamide's hallmark is not on carbon: one occurrence, no key atoms.
    EXPECT_EQ(key_atoms(compound_smiles("Methanesulfonamide"), Group::Sulfonamide), (Keys{{}}));
    EXPECT_TRUE(arise::chem::hallmarks(Group::Sulfonamide).empty());
    EXPECT_EQ(arise::chem::hallmarks(Group::CarboxylicEster).size(), 2U);
}

TEST(Groups, NamesAreReadable) {
    std::set<std::string_view> names;
    for (const auto group : arise::chem::kAllGroups) {
        const auto name = arise::chem::to_string(group);
        EXPECT_FALSE(name.empty());
        EXPECT_TRUE(names.insert(name).second) << "duplicate name " << name;
    }
    EXPECT_EQ(arise::chem::to_string(Group::CarboxylicAcid), "carboxylic acid");
    EXPECT_EQ(arise::chem::to_string(Group::Methoxy), "methoxy");
    EXPECT_EQ(arise::chem::kAllGroups.size(), 28U);  // the signed-off definitions table
}

// ---- exchangeable hydrogens --------------------------------------------------------------
// H on O, N, S and Se (the chemist's definition, 2026-11-07): the protons a D2O shake
// removes. Each SMILES's per-atom H counts are also checked against RDKit by the spec
// oracle (tests/data/reference/spec_rdkit_oracle.tsv).

int exchangeable(std::string_view smiles) {
    return arise::chem::exchangeable_hydrogens(parse_ok(smiles));
}

TEST(ExchangeableHydrogens, CountsHydrogenOnOxygen) {
    EXPECT_EQ(exchangeable(compound_smiles("Ethanol")), 1);          // ethanol
    EXPECT_EQ(exchangeable(compound_smiles("Acetic Acid")), 1);      // acetic acid
    EXPECT_EQ(exchangeable(compound_smiles("Ethylene Glycol")), 2);  // ethylene glycol
    EXPECT_EQ(exchangeable("O"), 2);                                 // water
}

TEST(ExchangeableHydrogens, CountsHydrogenOnNitrogenSulfurAndSelenium) {
    EXPECT_EQ(exchangeable(compound_smiles("Glycine")), 3);        // glycine: NH2 and OH
    EXPECT_EQ(exchangeable(compound_smiles("DL-Serine")), 4);      // serine: NH2, CH2-OH and COOH
    EXPECT_EQ(exchangeable(compound_smiles("Dimethylamine")), 1);  // dimethylamine
    EXPECT_EQ(exchangeable(compound_smiles("Formamide")), 2);      // formamide
    EXPECT_EQ(exchangeable(compound_smiles("Methanethiol")), 1);   // methanethiol
    EXPECT_EQ(exchangeable(compound_smiles("Methylselenol")), 1);  // methaneselenol
}

TEST(ExchangeableHydrogens, IgnoresHydrogenOnCarbonPhosphorusSiliconAndBoron) {
    EXPECT_EQ(exchangeable(compound_smiles("Dimethyl Ether")), 0);  // dimethyl ether
    EXPECT_EQ(exchangeable(compound_smiles("Trimethylamine")), 0);  // trimethylamine
    EXPECT_EQ(exchangeable(compound_smiles("Methylphosphine")),
              0);  // methylphosphine: P-H does not exchange
    EXPECT_EQ(exchangeable(compound_smiles("Silaethane")), 0);    // methylsilane
    EXPECT_EQ(exchangeable(compound_smiles("methylborane")), 0);  // methylborane
}

TEST(ExchangeableHydrogens, CountsExplicitHydrogenAtomsOnce) {
    // An explicit [H] atom bonded to O is that O's hydrogen, not a second one.
    EXPECT_EQ(exchangeable(syntax_case("explicit hydrogen atom on methanol O")), 1);
    EXPECT_EQ(exchangeable(syntax_case("explicit hydrogen atoms on methanol C")), 1);
}

}  // namespace
