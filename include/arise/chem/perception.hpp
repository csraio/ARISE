#pragma once

#include "arise/chem/molecule.hpp"
#include "arise/nmr1d/dept.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

namespace arise::chem {

/// The DEPT type of carbon atom `atom` (C, CH, CH2, CH3, by attached hydrogens);
/// empty for non-carbon atoms and for a carbon with four hydrogens (methane).
[[nodiscard]] std::optional<nmr1d::CarbonType> carbon_type(const Molecule& molecule,
                                                           std::size_t atom) noexcept;

/// Six-membered rings of C or N whose ring bonds alternate single/double: Kekulé benzene
/// and pyridine rings, also within fused systems. Each ring lists its atoms in ascending
/// order; rings are sorted, so the result does not depend on how the SMILES was written.
[[nodiscard]] std::vector<std::array<std::size_t, 6>> aromatic_rings(const Molecule& molecule);

/// Topological equivalence class of every atom: equal numbers mean equivalent atoms
/// (for carbons, the same expected 13C environment). Classes are numbered from 0 in
/// order of first appearance. Computed by iterative refinement of atom invariants
/// (element, isotope, charge, total hydrogens, heavy neighbours) with bonds of aromatic
/// rings treated as equivalent, so Kekulé benzene stays symmetric. Explicit [H] atoms count
/// only toward their neighbour's hydrogens, so isotope labels on H (D) are not
/// distinguished. Refinement finds true equivalence
/// in ordinary organic molecules; in rare highly regular graphs it can merge atoms that
/// are not equivalent.
[[nodiscard]] std::vector<std::size_t> symmetry_classes(const Molecule& molecule);

/// Functional groups a chemist can confirm or rule out.
enum class Group : std::uint8_t {
    CarboxylicAcid,  // C(=O)OH
    Aldehyde,        // C(=O)H, the carbonyl carbon bonded to no O, N or S
    Ketone,          // C(=O) bonded to exactly two carbons
    Ester,           // C(=O)O-C, including lactones
    Amide,           // C(=O)N
    Alcohol,         // OH on a carbon with only single bonds
    Phenol,          // OH on an aromatic-ring carbon
    Ether,           // C-O-C, neither carbon a carbonyl carbon
    Amine,           // N with only single bonds, not bonded to a carbonyl carbon
    Nitrile,         // C#N
    Alkene,          // C=C outside aromatic rings
    BenzeneRing,     // an aromatic ring of six carbons
};

/// All groups, in enum order.
inline constexpr std::array kAllGroups{
    Group::CarboxylicAcid, Group::Aldehyde, Group::Ketone, Group::Ester,
    Group::Amide,          Group::Alcohol,  Group::Phenol, Group::Ether,
    Group::Amine,          Group::Nitrile,  Group::Alkene, Group::BenzeneRing,
};

/// "carboxylic acid", "aldehyde", ...
[[nodiscard]] std::string_view to_string(Group group) noexcept;

/// How many times `group` occurs in the molecule (per carbonyl carbon, per OH, per C=C,
/// per ring, ...).
[[nodiscard]] int count_group(const Molecule& molecule, Group group);

}  // namespace arise::chem
