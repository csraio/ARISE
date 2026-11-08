// ARISE's molecule reader against BMRB's experimental entries: every case's structure, as
// the entry's own atom and bond list gives it (tests/data/reference/bmrb_cases.tsv, built
// by scripts/build_bmrb_reference.py and checked there against the entry's InChI), must
// read in ARISE with each atom carrying the hydrogens the entry bonds to it.

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "bmrb_cases.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <string>

namespace {

using arise::test::bmrb_cases;
using arise::test::BmrbCase;

// "" when ARISE reads the case's SMILES with the entry's hydrogens, else what differs.
std::string structure_difference(const BmrbCase& bmrb) {
    const auto molecule = arise::chem::parse_smiles(bmrb.smiles);
    if (!molecule) {
        return "does not parse: " + molecule.error().message;
    }
    const auto& parsed = molecule.value();
    if (parsed.atoms().size() != bmrb.hydrogens.size()) {
        return std::to_string(parsed.atoms().size()) + " atoms, the entry has " +
               std::to_string(bmrb.hydrogens.size());
    }
    std::string differences;
    for (std::size_t atom = 0; atom < bmrb.hydrogens.size(); ++atom) {
        if (parsed.hydrogen_count(atom) != bmrb.hydrogens[atom]) {
            differences += " " + bmrb.atom_names[atom] + " has " +
                           std::to_string(parsed.hydrogen_count(atom)) + " H, the entry " +
                           std::to_string(bmrb.hydrogens[atom]);
        }
    }
    return differences;
}

TEST(BmrbReference, EveryCaseStructureReadsWithTheEntrysHydrogens) {
    ASSERT_FALSE(bmrb_cases().empty());
    for (const auto& bmrb : bmrb_cases()) {
        EXPECT_EQ(structure_difference(bmrb), "") << bmrb.id() << " " << bmrb.smiles;
    }
}

// "" when every 13C-dimension assignment names a carbon and every 1H-dimension one an atom
// that carries H, in ARISE's reading of the structure; else the assignments that do not.
std::string assignment_difference(const BmrbCase& bmrb) {
    const auto molecule = arise::chem::parse_smiles(bmrb.smiles);
    if (!molecule) {
        return "does not parse";
    }
    std::string differences;
    for (const auto& peak : bmrb.peaks) {
        for (const auto atom : peak.c_atoms) {
            if (molecule.value().atoms()[atom].element != arise::chem::Element::C) {
                differences += " " + peak.list + " 13C on " + bmrb.atom_names[atom];
            }
        }
        for (const auto atom : peak.h_atoms) {
            if (molecule.value().hydrogen_count(atom) == 0) {
                differences += " " + peak.list + " 1H on " + bmrb.atom_names[atom];
            }
        }
    }
    return differences;
}

TEST(BmrbReference, AssignedAtomsAreTheListsNucleus) {
    for (const auto& bmrb : bmrb_cases()) {
        EXPECT_EQ(assignment_difference(bmrb), "") << bmrb.id();
    }
}

TEST(BmrbReference, ExchangeableHydrogensAreTheEntrysHOnONSAndSe) {
    // The entry's own bonds put each H on an atom; chem::exchangeable_hydrogens must count
    // exactly those on O, N, S and Se.
    for (const auto& bmrb : bmrb_cases()) {
        const auto molecule = arise::chem::parse_smiles(bmrb.smiles);
        ASSERT_TRUE(molecule) << bmrb.id();
        int expected = 0;
        for (std::size_t atom = 0; atom < bmrb.hydrogens.size(); ++atom) {
            const auto element = molecule.value().atoms()[atom].element;
            if (element == arise::chem::Element::O || element == arise::chem::Element::N ||
                element == arise::chem::Element::S || element == arise::chem::Element::Se) {
                expected += bmrb.hydrogens[atom];
            }
        }
        EXPECT_EQ(arise::chem::exchangeable_hydrogens(molecule.value()), expected) << bmrb.id();
    }
}

}  // namespace
