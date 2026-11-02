#include "arise/chem/perception.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/io/pubchem.hpp"
#include "arise/nmr1d/dept.hpp"
#include "pubchem_fixture.hpp"

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

using arise::chem::aromatic_rings;
using arise::chem::carbon_type;
using arise::chem::count_group;
using arise::chem::Element;
using arise::chem::Group;
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
    const auto molecule = parse_ok("CCOC(=O)C");  // ethyl acetate
    EXPECT_EQ(carbon_type(molecule, 0), CarbonType::CH3);
    EXPECT_EQ(carbon_type(molecule, 1), CarbonType::CH2);
    EXPECT_EQ(carbon_type(molecule, 2), std::nullopt);  // oxygen
    EXPECT_EQ(carbon_type(molecule, 3), CarbonType::C);
    EXPECT_EQ(carbon_type(molecule, 5), CarbonType::CH3);
    EXPECT_EQ(carbon_type(parse_ok("C1=CC=CC=C1"), 0), CarbonType::CH);
}

TEST(CarbonType, MethaneHasNoDeptType) {
    EXPECT_EQ(carbon_type(parse_ok("C"), 0), std::nullopt);
}

TEST(CarbonType, ExplicitHydrogenAtomsCount) {
    EXPECT_EQ(carbon_type(parse_ok("[H]C([H])(C)C"), 1), CarbonType::CH2);
}

// ---- Aromatic rings --------------------------------------------------------------

TEST(AromaticRings, FindsKekuleBenzene) {
    const auto rings = aromatic_rings(parse_ok("CC1=CC=CC=C1"));  // toluene
    ASSERT_EQ(rings.size(), 1U);
    EXPECT_EQ(rings[0], (std::vector<std::size_t>{1, 2, 3, 4, 5, 6}));
}

TEST(AromaticRings, FindsPyridine) {
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CC=NC=C1")).size(), 1U);
}

TEST(AromaticRings, FindsBothRingsOfKekuleNaphthalene) {
    // Written this way, the right-hand ring has only two double bonds of its own; the
    // third is the fusion bond's partner in the left ring. Both rings are aromatic.
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CC=C2C=CC=CC2=C1")).size(), 2U);
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CC=C(C=C1)C2=CC=CC=C2")).size(), 2U);  // biphenyl
}

TEST(AromaticRings, FindsFiveMemberedHeteroaromatics) {
    // Huckel 6 pi: two ring C=C plus the lone pair of O, S or N(H/R) (user decision:
    // furan, pyrrole and thiophene are aromatic).
    EXPECT_EQ(aromatic_rings(parse_ok("C1=COC=C1")).size(), 1U);   // furan
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CNC=C1")).size(), 1U);   // pyrrole
    EXPECT_EQ(aromatic_rings(parse_ok("CN1C=CC=C1")).size(), 1U);  // N-methylpyrrole
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CSC=C1")).size(), 1U);   // thiophene
    EXPECT_EQ(aromatic_rings(parse_ok("C1=COC=C1"))[0].size(), 5U);
    // Indole: the benzene ring and the pyrrole ring, which shares a fusion bond. In both
    // Kekule forms the fusion carbons' double bonds lie in the benzene ring.
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CC=C2C(=C1)C=CN2")).size(), 2U);
    EXPECT_EQ(aromatic_rings(parse_ok("N1C=CC2=CC=CC=C12")).size(), 2U);
}

TEST(AromaticRings, IgnoresNonAromaticFiveRings) {
    EXPECT_TRUE(aromatic_rings(parse_ok("C1COC=C1")).empty());   // 2,3-dihydrofuran
    EXPECT_TRUE(aromatic_rings(parse_ok("C1C=CC=C1")).empty());  // cyclopentadiene: no lone pair
    EXPECT_TRUE(aromatic_rings(parse_ok("C1CCOC1")).empty());    // THF
    // A ring S with two extra bonds has no free lone pair for the ring.
    EXPECT_TRUE(aromatic_rings(parse_ok("CS1(C)C=CC=C1")).empty());
    // Fluoranthene: three benzenoid rings; its five-membered ring has five carbons and no
    // lone-pair donor, so it is not counted.
    EXPECT_EQ(aromatic_rings(parse_ok("C1=CC2=CC=CC3=C2C(=C1)C4=CC=CC=C43")).size(), 3U);
}

TEST(AromaticRings, IgnoresNonAromaticSixRings) {
    EXPECT_TRUE(aromatic_rings(parse_ok("C1=CCCCC1")).empty());         // cyclohexene
    EXPECT_TRUE(aromatic_rings(parse_ok("C1=CC=CCC1")).empty());        // cyclohexadiene
    EXPECT_TRUE(aromatic_rings(parse_ok("C1CCCCC1")).empty());          // cyclohexane
    EXPECT_TRUE(aromatic_rings(parse_ok("C1=CC(=O)C=CC1=O")).empty());  // benzoquinone
    EXPECT_TRUE(aromatic_rings(parse_ok("CS1=CC=CC=C1")).empty());      // thiabenzene ylide
}

TEST(AromaticRings, DoNotDependOnHowTheSmilesIsWritten) {
    const auto a = aromatic_rings(parse_ok("C1=CC=C(C=C1)C"));
    const auto b = aromatic_rings(parse_ok("C1=CC(=CC=C1)C"));
    ASSERT_EQ(a.size(), 1U);
    ASSERT_EQ(b.size(), 1U);
    EXPECT_EQ(a[0], b[0]);  // both rings are atoms 0-5
}

// ---- Symmetry classes ------------------------------------------------------------

TEST(Symmetry, NumbersClassesInOrderOfFirstAppearance) {
    // Toluene: CH3, ipso, ortho, meta, para, meta, ortho.
    EXPECT_EQ(symmetry_classes(parse_ok("CC1=CC=CC=C1")),
              (std::vector<std::size_t>{0, 1, 2, 3, 4, 3, 2}));
}

TEST(Symmetry, KekuleBenzeneStaysSymmetric) {
    EXPECT_EQ(carbon_class_sizes("C1=CC=CC=C1"), (std::vector<int>{6}));
    // o-Xylene: the two ortho CH and the two ipso C must pair up although one ipso
    // carbon has a written double bond to its neighbour and the other a single bond.
    EXPECT_EQ(carbon_class_sizes("CC1=CC=CC=C1C"), (std::vector<int>{2, 2, 2, 2}));
}

TEST(Symmetry, RegioisomersGiveTheirTextbookCarbonCounts) {
    // Distinct 13C environments: toluene 5, o-/m-/p-xylene 4/5/3, ethylbenzene 6.
    EXPECT_EQ(carbon_class_sizes("CC1=CC=CC=C1").size(), 5U);
    EXPECT_EQ(carbon_class_sizes("CC1=CC=CC=C1C").size(), 4U);
    EXPECT_EQ(carbon_class_sizes("CC1=CC(C)=CC=C1"), (std::vector<int>{1, 1, 2, 2, 2}));
    EXPECT_EQ(carbon_class_sizes("CC1=CC=C(C)C=C1"), (std::vector<int>{2, 2, 4}));
    EXPECT_EQ(carbon_class_sizes("CCC1=CC=CC=C1").size(), 6U);
}

TEST(Symmetry, AliphaticEquivalences) {
    EXPECT_EQ(carbon_class_sizes("CC(C)O"), (std::vector<int>{1, 2}));     // isopropanol
    EXPECT_EQ(carbon_class_sizes("CCOCC"), (std::vector<int>{2, 2}));      // diethyl ether
    EXPECT_EQ(carbon_class_sizes("CC(C)(C)C"), (std::vector<int>{1, 4}));  // neopentane
    EXPECT_EQ(carbon_class_sizes("CCOC(=O)C"), (std::vector<int>{1, 1, 1, 1}));
}

TEST(Symmetry, FiveMemberedHeteroaromaticsStaySymmetric) {
    EXPECT_EQ(carbon_class_sizes("C1=COC=C1"), (std::vector<int>{2, 2}));  // furan
    EXPECT_EQ(carbon_class_sizes("C1=CSC=C1"), (std::vector<int>{2, 2}));  // thiophene
    EXPECT_EQ(carbon_class_sizes("CC1=CC=C(C)O1"),
              (std::vector<int>{2, 2, 2}));  // 2,5-dimethylfuran
}

TEST(Symmetry, NaphthaleneHasThreeCarbonEnvironments) {
    // alpha (4), beta (4) and the two fusion carbons.
    EXPECT_EQ(carbon_class_sizes("C1=CC=C2C=CC=CC2=C1"), (std::vector<int>{2, 4, 4}));
}

TEST(Symmetry, DoesNotDependOnHowTheSmilesIsWritten) {
    EXPECT_EQ(carbon_class_sizes("CC1=CC=CC=C1"), carbon_class_sizes("C1=CC=C(C=C1)C"));
    EXPECT_EQ(carbon_class_sizes("CC1=CC=C(C)C=C1"), carbon_class_sizes("C1=CC(=CC=C1C)C"));
}

TEST(Symmetry, ExplicitHydrogensDoNotChangeEquivalence) {
    // Ethane with one methyl's hydrogens written as atoms: both CH3 are still equivalent.
    EXPECT_EQ(carbon_class_sizes("[H]C([H])([H])C"), (std::vector<int>{2}));
}

TEST(Symmetry, HydrogenCountDistinguishesOtherwiseAlikeAtoms) {
    // Propane with a CH2 radical end: the two end carbons differ only in hydrogens.
    EXPECT_EQ(carbon_class_sizes("[CH2]CC"), (std::vector<int>{1, 1, 1}));
}

TEST(Symmetry, ResonanceSharedOxygensAreEquivalent) {
    // The two O of a nitro group or a carboxylate share the charge by resonance.
    const auto nitro = parse_ok("CC[N+](=O)[O-]");
    const auto classes = arise::chem::symmetry_classes(nitro);
    EXPECT_EQ(classes[3], classes[4]);
    const auto glycine = parse_ok("C(C(=O)[O-])[NH3+]");
    const auto zwitterion = arise::chem::symmetry_classes(glycine);
    EXPECT_EQ(zwitterion[2], zwitterion[3]);
    // ...but an acid's OH and C=O stay distinct.
    const auto acid = arise::chem::symmetry_classes(parse_ok("CC(=O)O"));
    EXPECT_NE(acid[2], acid[3]);
}

TEST(Symmetry, IsotopeLabelBreaksEquivalence) {
    // 13C-labelled methyl in propane: the two methyls are no longer equivalent.
    EXPECT_EQ(carbon_class_sizes("[13CH3]CC"), (std::vector<int>{1, 1, 1}));
}

// ---- NMR equivalence ---------------------------------------------------------------
// Expected 13C environment counts. Literature behaviour (e.g. Pretsch et al.; standard
// texts on diastereotopicity and amide rotation): valine's isopropyl methyls differ,
// linalool's =C(CH3)2 methyls differ (~17.6/25.7 ppm), DMF's N-methyls differ.

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
    EXPECT_EQ(topological("C=CC(C)(O)CCC=C(C)C"), 9);  // linalool
    EXPECT_EQ(environments("C=CC(C)(O)CCC=C(C)C"), (Range{10, 10}));
    EXPECT_EQ(environments("CC=C(C)C"), (Range{5, 5}));     // 2-methyl-2-butene
    EXPECT_EQ(environments("ON=C1CCCCC1"), (Range{6, 6}));  // cyclohexanone oxime
    EXPECT_EQ(environments("CC(C)=NO"), (Range{3, 3}));     // acetone oxime
}

TEST(NmrEquivalence, EZNeedsDifferentGroupsOnTheOtherEnd) {
    EXPECT_EQ(environments("CC(C)=C(C)C"), (Range{2, 2}));  // tetramethylethylene
    EXPECT_EQ(environments("C=C(C)C"), (Range{3, 3}));      // isobutylene: =CH2 end
}

TEST(NmrEquivalence, AmideRotationSplitsAlways) {
    EXPECT_EQ(topological("CN(C)C=O"), 2);  // DMF
    EXPECT_EQ(environments("CN(C)C=O"), (Range{3, 3}));
    EXPECT_EQ(environments("CC(=O)N(C)C"), (Range{4, 4}));  // N,N-dimethylacetamide
    EXPECT_EQ(environments("CN(C)C=S"), (Range{3, 3}));     // thioamide
}

TEST(NmrEquivalence, BorderlineRotationGivesARange) {
    EXPECT_EQ(environments("CN(C)C(=O)N(C)C"), (Range{2, 3}));  // tetramethylurea
    EXPECT_EQ(environments("COC(=O)N(C)C"), (Range{3, 4}));     // methyl dimethylcarbamate
    EXPECT_EQ(environments("CN(C)C=CC(C)=O"), (Range{5, 6}));   // enaminone
}

TEST(NmrEquivalence, DiastereotopicGroupsGiveARange) {
    EXPECT_EQ(environments("CC(C)C(N)C(=O)O"), (Range{4, 5}));           // valine
    EXPECT_EQ(environments("CC(C)CC(C)O"), (Range{5, 6}));               // 4-methyl-2-pentanol
    EXPECT_EQ(environments("C=CS(=O)C(C)C"), (Range{4, 5}));             // chiral sulfoxide
    EXPECT_EQ(environments("CP(=O)(C(C)C)C1=CC=CC=C1"), (Range{7, 8}));  // chiral P
}

TEST(NmrEquivalence, SeveralProchiralCentresAlsoGiveARange) {
    // nmrshiftdb2 10017805: no stereocentre, yet each acetal's two OEt groups give
    // different peaks, because several prochiral centres make each other's groups
    // diastereotopic. Over-splitting only widens the range, which is the safe direction.
    const auto range = environments("CCOC(OCC)C(C(C)C)C(OCC)OCC");
    EXPECT_EQ(range.fewest, 6);
    EXPECT_GE(range.most, 8);
}

TEST(NmrEquivalence, OneProchiralCentreAloneKeepsEquivalence) {
    // Isopropanol: one prochiral centre, no stereocentre; the methyls are enantiotopic
    // and give one peak.
    EXPECT_EQ(environments("CC(C)O"), (Range{2, 2}));
}

TEST(NmrEquivalence, UnknownStereoCanBreakASymmetry) {
    // Hexitols: mannitol shows 3 carbon environments, sorbitol 6; SMILES without stereo
    // cannot tell them apart.
    EXPECT_EQ(environments("OCC(O)C(O)C(O)C(O)CO"), (Range{3, 6}));
}

TEST(NmrEquivalence, FastRotationKeepsEquivalence) {
    // tert-butyl methyls (three identical groups) and a phenyl's ortho carbons stay
    // equivalent even next to a stereocentre.
    EXPECT_EQ(environments("CC(C)(C)C(C)O"), (Range{4, 4}));     // 3,3-dimethyl-2-butanol
    EXPECT_EQ(environments("CC(O)C1=CC=CC=C1"), (Range{6, 6}));  // 1-phenylethanol
    EXPECT_EQ(environments("CCOC(=O)C"), (Range{4, 4}));         // ethyl acetate
    EXPECT_EQ(environments("CC1=CC=C(C)C=C1"), (Range{3, 3}));   // p-xylene
}

TEST(NmrEquivalence, ProchiralCentresBeyondCarbonMaySplit) {
    // A quaternary N+ does not invert, so with a stereocentre elsewhere its two methyls
    // are diastereotopic (the "most" partition splits them): 6..7 carbon environments.
    EXPECT_EQ(environments("C[N+](C)(CC)CC([O-])C"), (Range{6, 7}));
    // A neutral amine inverts fast, so its two methyls stay equivalent.
    EXPECT_EQ(environments("CN(C)CC(O)C"), (Range{4, 4}));
}

TEST(NmrEquivalence, NeverFewerClassesThanTopology) {
    for (const auto* smiles : {"CC(C)C(N)C(=O)O", "CN(C)C=O", "OCC(O)C(O)C(O)C(O)CO"}) {
        const auto range = environments(smiles);
        EXPECT_GE(range.fewest, topological(smiles)) << smiles;
        EXPECT_GE(range.most, range.fewest) << smiles;
    }
}

// ---- Functional groups (docs/chem_group_definitions.md) -----------------------------
// One test per family; each checks the inclusions and exclusions the definitions table
// names. External checks against ChEBI and RDKit live in reference_test.cpp.

TEST(Groups, AcidsAldehydesKetones) {
    EXPECT_EQ(groups("CC(=O)O", Group::CarboxylicAcid), 1);
    EXPECT_EQ(groups("C(=O)O", Group::CarboxylicAcid), 1);             // formic acid...
    EXPECT_EQ(groups("C(=O)O", Group::Aldehyde), 0);                   // ...not an aldehyde
    EXPECT_EQ(groups("OC(=O)C1=CC=CC=C1", Group::CarboxylicAcid), 1);  // benzoic acid
    EXPECT_EQ(groups("CC=O", Group::Aldehyde), 1);
    EXPECT_EQ(groups("C=O", Group::Aldehyde), 1);       // formaldehyde (user decision)
    EXPECT_EQ(groups("COC=O", Group::Aldehyde), 0);     // methyl formate
    EXPECT_EQ(groups("CN(C)C=O", Group::Aldehyde), 0);  // DMF
    EXPECT_EQ(groups("CC(=O)C", Group::Ketone), 1);
    EXPECT_EQ(groups("CC(=O)C1=CC=CC=C1", Group::Ketone), 1);  // acetophenone
    EXPECT_EQ(groups("C1=CC(=O)C=CC1=O", Group::Ketone), 2);   // benzoquinone
    EXPECT_EQ(groups("CCOC(=O)C", Group::Ketone), 0);
    EXPECT_EQ(groups("CC=O", Group::Ketone), 0);  // an aldehyde has one carbon side
}

TEST(Groups, EstersAndTheirNeighbours) {
    EXPECT_EQ(groups("CCOC(=O)C", Group::CarboxylicEster), 1);      // ethyl acetate
    EXPECT_EQ(groups("COC=O", Group::CarboxylicEster), 1);          // formate (user decision)
    EXPECT_EQ(groups("C1CC(=O)OC1", Group::CarboxylicEster), 1);    // lactone
    EXPECT_EQ(groups("COC(=O)OC", Group::CarboxylicEster), 0);      // carbonate
    EXPECT_EQ(groups("COC(N)=O", Group::CarboxylicEster), 0);       // carbamate
    EXPECT_EQ(groups("CC(=O)OC(C)=O", Group::CarboxylicEster), 0);  // anhydride
    EXPECT_EQ(groups("COC(=O)OC", Group::CarbonateEster), 1);
    EXPECT_EQ(groups("O=C1OCCO1", Group::CarbonateEster), 1);  // ethylene carbonate
    EXPECT_EQ(groups("CCOC(=O)C", Group::CarbonateEster), 0);
    EXPECT_EQ(groups("COC(=O)O", Group::CarbonateEster), 0);  // monoester: one O is an OH
    EXPECT_EQ(groups("CC(=O)OC(C)=O", Group::CarboxylicAnhydride), 1);
    EXPECT_EQ(groups("O=C1CCC(=O)O1", Group::CarboxylicAnhydride), 1);  // succinic anhydride
    EXPECT_EQ(groups("CCOC(=O)C", Group::CarboxylicAnhydride), 0);
}

TEST(Groups, AmidesCarbamatesUreas) {
    EXPECT_EQ(groups("CC(=O)N", Group::Carboxamide), 1);
    EXPECT_EQ(groups("CN(C)C=O", Group::Carboxamide), 1);  // formamide
    EXPECT_EQ(groups("CNC(C)=O", Group::Carboxamide), 1);
    EXPECT_EQ(groups("O=C1CCCN1", Group::Carboxamide), 1);  // lactam
    EXPECT_EQ(groups("NC(N)=O", Group::Carboxamide), 0);    // urea (user decision)
    EXPECT_EQ(groups("COC(N)=O", Group::Carboxamide), 0);   // carbamate (user decision)
    EXPECT_EQ(groups("COC(N)=O", Group::CarbamateEster), 1);
    EXPECT_EQ(groups("O=C1NCCO1", Group::CarbamateEster), 1);  // oxazolidinone
    EXPECT_EQ(groups("NC(N)=O", Group::CarbamateEster), 0);
    EXPECT_EQ(groups("NC(=O)O", Group::CarbamateEster), 0);  // carbamic acid: no ester O
    // O-carbamoyl oximes count (user decision 2026-10-31, as in ChEBI): methomyl, aldicarb.
    EXPECT_EQ(groups("CNC(=O)ON=C(C)SC", Group::CarbamateEster), 1);
    EXPECT_EQ(groups("CNC(=O)ON=CC(C)(C)SC", Group::CarbamateEster), 1);
    EXPECT_EQ(groups("CNC(=O)ON=C(C)SC", Group::CarboxylicEster), 0);
    EXPECT_EQ(groups("CNC(=O)ONC", Group::CarbamateEster), 0);       // O on an amine N: not oxime
    EXPECT_EQ(groups("CC(=O)ON=C(C)C", Group::CarboxylicEster), 0);  // O-acyl oxime: not decided
    EXPECT_EQ(groups("NC(N)=O", Group::Urea), 1);
    EXPECT_EQ(groups("CN(C)C(=O)N(C)C", Group::Urea), 1);
    EXPECT_EQ(groups("O=C1NCCN1", Group::Urea), 1);  // imidazolidinone
    EXPECT_EQ(groups("CC(=O)N", Group::Urea), 0);
}

TEST(Groups, HydroxylFamilies) {
    EXPECT_EQ(groups("CCO", Group::Alcohol), 1);
    EXPECT_EQ(groups("C=CO", Group::Alcohol), 0);
    EXPECT_EQ(groups("OC1=CC=CC=C1", Group::Alcohol), 0);
    EXPECT_EQ(groups("CC(=O)O", Group::Alcohol), 0);
    EXPECT_EQ(groups("C=CO", Group::Enol), 1);
    EXPECT_EQ(groups("CCO", Group::Enol), 0);
    EXPECT_EQ(groups("OC1=CC=CC=C1", Group::Enol), 0);
    EXPECT_EQ(groups("CC=C=O", Group::Enol), 0);  // ketene: the O carries no H
    EXPECT_EQ(groups("OC1=CC=CC=C1", Group::Phenol), 1);
    EXPECT_EQ(groups("OC1=CC2=CC=CC=C2C=C1", Group::Phenol), 1);  // 2-naphthol
    EXPECT_EQ(groups("OC(=O)C1=CC=C(O)C=C1", Group::Phenol), 1);
    EXPECT_EQ(groups("OC1=CN=CC=C1", Group::Phenol), 0);  // hydroxypyridine
    EXPECT_EQ(groups("OC1=CN=CC=C1", Group::Pyridone), 1);
    EXPECT_EQ(groups("O=C1C=CC=CN1", Group::Pyridone), 1);  // 2-pyridone
    EXPECT_EQ(groups("O=C1C=CNC=C1", Group::Pyridone), 1);  // 4-pyridone
    EXPECT_EQ(groups("OC1=CC=CC=C1", Group::Pyridone), 0);
    EXPECT_EQ(groups("O=C1CCCCN1", Group::Pyridone), 0);  // valerolactam: no ring C=C
}

TEST(Groups, Ethers) {
    EXPECT_EQ(groups("CCOCC", Group::Ether), 1);
    EXPECT_EQ(groups("COCC(=O)C", Group::Ether), 1);      // methoxyacetone
    EXPECT_EQ(groups("COC1=CC=CC=C1", Group::Ether), 1);  // anisole
    EXPECT_EQ(groups("C1COCCO1", Group::Ether), 2);       // dioxane
    EXPECT_EQ(groups("C1CCOC1", Group::Ether), 1);        // THF
    EXPECT_EQ(groups("CC1CO1", Group::Ether), 1);         // epoxide
    EXPECT_EQ(groups("CCOC(=O)C", Group::Ether), 0);      // ester O
    EXPECT_EQ(groups("C1=COC=C1", Group::Ether), 0);      // furan O (user decision)
}

TEST(Groups, NitrogenFamilies) {
    EXPECT_EQ(groups("CN", Group::Amine), 1);
    EXPECT_EQ(groups("CCN(CC)CC", Group::Amine), 1);
    EXPECT_EQ(groups("NC1=CC=CC=C1", Group::Amine), 1);  // aniline
    EXPECT_EQ(groups("CC(=O)N", Group::Amine), 0);
    EXPECT_EQ(groups("C1=CC=NC=C1", Group::Amine), 0);  // pyridine N
    EXPECT_EQ(groups("C1=CNC=C1", Group::Amine), 0);    // pyrrole N (user decision)
    EXPECT_EQ(groups("CS(N)(=O)=O", Group::Amine), 0);  // sulfonamide N
    EXPECT_EQ(groups("CNN", Group::Amine), 0);          // hydrazine N
    EXPECT_EQ(groups("CNO", Group::Amine), 0);          // hydroxylamine N
    EXPECT_EQ(groups("CC(N)=S", Group::Amine), 0);      // thioamide N
    EXPECT_EQ(groups("CS(N)(=O)=O", Group::Sulfonamide), 1);
    EXPECT_EQ(groups("CN(C)S(=O)(=O)C", Group::Sulfonamide), 1);
    EXPECT_EQ(groups("CS(N)=O", Group::Sulfonamide), 0);  // sulfinamide: one S=O
    EXPECT_EQ(groups("CNN", Group::Hydrazine), 1);
    EXPECT_EQ(groups("CC(=O)NN", Group::Hydrazine), 1);  // hydrazide
    EXPECT_EQ(groups("CNO", Group::Hydroxylamine), 1);
    EXPECT_EQ(groups("CC(C)=NO", Group::Hydroxylamine), 0);  // oxime
    EXPECT_EQ(groups("CC(C)=NC", Group::Imine), 1);
    EXPECT_EQ(groups("CC(C)=NO", Group::Imine), 0);
    EXPECT_EQ(groups("CC(C)=NN", Group::Imine), 0);  // hydrazone
    EXPECT_EQ(groups("CC(C)=NO", Group::Oxime), 1);
    EXPECT_EQ(groups("CC(C)=NC", Group::Oxime), 0);
    EXPECT_EQ(groups("CC(C)=NOC", Group::Oxime), 0);  // O-methyl oxime: no OH
    EXPECT_EQ(groups("CC#N", Group::Nitrile), 1);
}

TEST(Groups, AlkenesAndAromaticRings) {
    EXPECT_EQ(groups("C=CC=C", Group::Alkene), 2);
    EXPECT_EQ(groups("C=CC1=CC=CC=C1", Group::Alkene), 1);  // styrene: ring bonds excluded
    EXPECT_EQ(groups("C1=CC=CC=C1", Group::Alkene), 0);
    EXPECT_EQ(groups("C1=CC(=O)C=CC1=O", Group::Alkene), 2);  // quinone is not aromatic
    EXPECT_EQ(groups("C1=COC=C1", Group::Alkene), 0);         // furan ring bonds
    EXPECT_EQ(groups("C=CC1=CC=CO1", Group::Alkene), 1);      // 2-vinylfuran
    EXPECT_EQ(groups("C=CC1=CC=CC=C1", Group::BenzeneRing), 1);
    EXPECT_EQ(groups("C1=CC=C2C=CC=CC2=C1", Group::BenzeneRing), 2);
    EXPECT_EQ(groups("C1=CC=NC=C1", Group::BenzeneRing), 0);
    EXPECT_EQ(groups("C1=CC=NC=C1", Group::PyridineRing), 1);
    EXPECT_EQ(groups("C1=CC=C2C(=C1)C=CC=N2", Group::PyridineRing), 1);  // quinoline
    EXPECT_EQ(groups("C1=CC=C2C(=C1)C=CC=N2", Group::BenzeneRing), 1);
    EXPECT_EQ(groups("C1=COC=C1", Group::FuranRing), 1);
    EXPECT_EQ(groups("C1=CNC=C1", Group::FuranRing), 0);  // the heteroatom decides the ring
    EXPECT_EQ(groups("C1=CSC=C1", Group::PyrroleRing), 0);
    EXPECT_EQ(groups("C1=CC=C2C(=C1)C=CO2", Group::FuranRing), 1);  // benzofuran
    EXPECT_EQ(groups("C1CCOC1", Group::FuranRing), 0);              // THF
    EXPECT_EQ(groups("C1=CNC=C1", Group::PyrroleRing), 1);
    EXPECT_EQ(groups("CN1C=CC=C1", Group::PyrroleRing), 1);
    EXPECT_EQ(groups("C1=CC=C2C(=C1)C=CN2", Group::PyrroleRing), 1);  // indole
    EXPECT_EQ(groups("C1=CC=C2C(=C1)C=CN2", Group::BenzeneRing), 1);
    EXPECT_EQ(groups("C1CCNC1", Group::PyrroleRing), 0);  // pyrrolidine
    EXPECT_EQ(groups("C1=CSC=C1", Group::ThiopheneRing), 1);
    EXPECT_EQ(groups("C1CCSC1", Group::ThiopheneRing), 0);
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
            .groups = {{Group::CarboxylicEster, 1}},
        },
        {
            .cid = "CID 22172",
            .carbon_classes = {1, 1, 1, 1},
            .groups = {{Group::Ketone, 1}, {Group::Ether, 1}},
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
    EXPECT_EQ(groups("C(C(=O)[O-])[NH3+]", Group::CarboxylicAcid), 1);  // glycine zwitterion
    EXPECT_EQ(groups("C(C(=O)[O-])[NH3+]", Group::Amine), 1);
    EXPECT_EQ(groups("C(C[NH3+])C[O-]", Group::Alcohol), 1);
    EXPECT_EQ(groups("C[N+](C)(C)[O-]", Group::Amine), 0);  // trimethylamine N-oxide
    EXPECT_EQ(groups("C[N+](C)(C)[O-]", Group::Hydroxylamine), 0);
    EXPECT_EQ(groups("[C-]#[N+]C", Group::Nitrile), 0);  // methyl isocyanide
    EXPECT_EQ(groups("CC#N", Group::Nitrile), 1);
    EXPECT_EQ(groups("CCC(=[O+][O-])C", Group::Ketone), 0);             // carbonyl oxide
    EXPECT_EQ(groups("C1=CC=[N+](C=C1)[O-]", Group::PyridineRing), 1);  // pyridine N-oxide
}

TEST(Groups, NamesAreReadable) {
    std::set<std::string_view> names;
    for (const auto group : arise::chem::kAllGroups) {
        const auto name = arise::chem::to_string(group);
        EXPECT_FALSE(name.empty());
        EXPECT_TRUE(names.insert(name).second) << "duplicate name " << name;
    }
    EXPECT_EQ(arise::chem::to_string(Group::CarboxylicAcid), "carboxylic acid");
    EXPECT_EQ(arise::chem::kAllGroups.size(), 27U);  // the signed-off definitions table
}

}  // namespace
