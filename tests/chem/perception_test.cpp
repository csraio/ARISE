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
    EXPECT_EQ(rings[0], (std::array<std::size_t, 6>{1, 2, 3, 4, 5, 6}));
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

TEST(AromaticRings, IgnoresNonAromaticSixRings) {
    EXPECT_TRUE(aromatic_rings(parse_ok("C1=CCCCC1")).empty());         // cyclohexene
    EXPECT_TRUE(aromatic_rings(parse_ok("C1=CC=CCC1")).empty());        // cyclohexadiene
    EXPECT_TRUE(aromatic_rings(parse_ok("C1CCCCC1")).empty());          // cyclohexane
    EXPECT_TRUE(aromatic_rings(parse_ok("C1=CC(=O)C=CC1=O")).empty());  // benzoquinone
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

TEST(NmrEquivalence, NeverFewerClassesThanTopology) {
    for (const auto* smiles : {"CC(C)C(N)C(=O)O", "CN(C)C=O", "OCC(O)C(O)C(O)C(O)CO"}) {
        const auto range = environments(smiles);
        EXPECT_GE(range.fewest, topological(smiles)) << smiles;
        EXPECT_GE(range.most, range.fewest) << smiles;
    }
}

// ---- Functional groups -------------------------------------------------------------

TEST(Groups, CarbonylFamilies) {
    EXPECT_EQ(groups("CC(=O)O", Group::CarboxylicAcid), 1);  // acetic acid
    EXPECT_EQ(groups("C(=O)O", Group::CarboxylicAcid), 1);   // formic acid...
    EXPECT_EQ(groups("C(=O)O", Group::Aldehyde), 0);         // ...is not an aldehyde
    EXPECT_EQ(groups("CC=O", Group::Aldehyde), 1);           // acetaldehyde
    EXPECT_EQ(groups("C=O", Group::Aldehyde), 1);            // formaldehyde
    EXPECT_EQ(groups("CC=O", Group::Ketone), 0);
    EXPECT_EQ(groups("CC(=O)C", Group::Ketone), 1);  // acetone
    EXPECT_EQ(groups("CC(=O)N", Group::Amide), 1);   // acetamide
    EXPECT_EQ(groups("CC(=O)N", Group::Amine), 0);
    // N-methylacetamide: its N has a carbonyl and a methyl neighbour; still no amine.
    EXPECT_EQ(groups("CNC(C)=O", Group::Amine), 0);
    EXPECT_EQ(groups("CNC(C)=O", Group::Amide), 1);
}

TEST(Groups, Esters) {
    EXPECT_EQ(groups("CCOC(=O)C", Group::Ester), 1);  // ethyl acetate
    EXPECT_EQ(groups("CCOC(=O)C", Group::Ether), 0);  // its O is not an ether O
    EXPECT_EQ(groups("CCOC(=O)C", Group::Ketone), 0);
    EXPECT_EQ(groups("COC=O", Group::Ester), 1);            // methyl formate...
    EXPECT_EQ(groups("COC=O", Group::Aldehyde), 0);         // ...is not an aldehyde
    EXPECT_EQ(groups("C1CC(=O)OC1", Group::Ester), 1);      // gamma-butyrolactone
    EXPECT_EQ(groups("CC(=O)OC(C)=O", Group::Ester), 0);    // acetic anhydride
    EXPECT_EQ(groups("CCOC(=O)CC(C)=O", Group::Ester), 1);  // ethyl acetoacetate
    EXPECT_EQ(groups("CCOC(=O)CC(C)=O", Group::Ketone), 1);
}

TEST(Groups, HydroxylAndEtherOxygens) {
    EXPECT_EQ(groups("CCO", Group::Alcohol), 1);
    EXPECT_EQ(groups("OC1=CC=CC=C1", Group::Phenol), 1);
    EXPECT_EQ(groups("OC1=CC=CC=C1", Group::Alcohol), 0);
    EXPECT_EQ(groups("C=CO", Group::Alcohol), 0);  // an enol is not an alcohol
    EXPECT_EQ(groups("CC(=O)O", Group::Alcohol), 0);
    EXPECT_EQ(groups("CCOCC", Group::Ether), 1);
    EXPECT_EQ(groups("COCC(=O)C", Group::Ether), 1);  // methoxyacetone
    EXPECT_EQ(groups("COCC(=O)C", Group::Ketone), 1);
    EXPECT_EQ(groups("COC1=CC=CC=C1", Group::Ether), 1);  // anisole
    EXPECT_EQ(groups("OC(=O)C1=CC=C(O)C=C1", Group::CarboxylicAcid), 1);
    EXPECT_EQ(groups("OC(=O)C1=CC=C(O)C=C1", Group::Phenol), 1);
}

TEST(Groups, NitrogenGroups) {
    EXPECT_EQ(groups("CN", Group::Amine), 1);
    EXPECT_EQ(groups("CCN(CC)CC", Group::Amine), 1);  // triethylamine
    EXPECT_EQ(groups("CC#N", Group::Nitrile), 1);
    EXPECT_EQ(groups("CC#N", Group::Amine), 0);
    EXPECT_EQ(groups("C1=CC=NC=C1", Group::Amine), 0);  // pyridine N is not an amine
}

TEST(Groups, AlkenesAndRings) {
    EXPECT_EQ(groups("C=CC=C", Group::Alkene), 2);
    EXPECT_EQ(groups("C=CC1=CC=CC=C1", Group::Alkene), 1);  // styrene: ring bonds excluded
    EXPECT_EQ(groups("C=CC1=CC=CC=C1", Group::BenzeneRing), 1);
    EXPECT_EQ(groups("C1=CC=CC=C1", Group::Alkene), 0);
    EXPECT_EQ(groups("C1=CC=NC=C1", Group::BenzeneRing), 0);  // pyridine has N
    EXPECT_EQ(groups("C1=CC=C2C=CC=CC2=C1", Group::BenzeneRing), 2);
    EXPECT_EQ(groups("C1=CC(=O)C=CC1=O", Group::Ketone), 2);  // benzoquinone
    EXPECT_EQ(groups("C1=CC(=O)C=CC1=O", Group::Alkene), 2);
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
        {.cid = "CID 8857", .carbon_classes = {1, 1, 1, 1}, .groups = {{Group::Ester, 1}}},
        {.cid = "CID 11124", .carbon_classes = {1, 1, 1, 1}, .groups = {{Group::Ester, 1}}},
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

TEST(Groups, NamesAreReadable) {
    std::set<std::string_view> names;
    for (const auto group : arise::chem::kAllGroups) {
        const auto name = arise::chem::to_string(group);
        EXPECT_FALSE(name.empty());
        EXPECT_TRUE(names.insert(name).second) << "duplicate name " << name;
    }
    EXPECT_EQ(arise::chem::to_string(Group::CarboxylicAcid), "carboxylic acid");
}

}  // namespace
