#pragma once

#include "arise/chem/molecule.hpp"
#include "arise/nmr1d/dept.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace arise::chem {

/// The DEPT type of carbon atom `atom` (C, CH, CH2, CH3, by attached hydrogens);
/// empty for non-carbon atoms and for a carbon with four hydrogens (methane).
[[nodiscard]] std::optional<nmr1d::CarbonType> carbon_type(const Molecule& molecule,
                                                           std::size_t atom) noexcept;

/// Hydrogens on O, N, S and Se (Molecule::hydrogen_count of each): the protons that
/// exchange with D2O (OH, NH, SH, SeH), so the number a chemist can confirm from a D2O shake
/// or from 1H peaks without an HSQC partner (hsqc::suggest_exchangeable_count). H on P, Si
/// or B does not exchange and is not counted. The chemist's definition, 2026-11-07.
[[nodiscard]] int exchangeable_hydrogens(const Molecule& molecule) noexcept;

/// Aromatic rings, read from Kekulé SMILES:
/// - six-membered rings of C/N where every ring atom has one double bond (benzene,
///   pyridine);
/// - five-membered rings of four C plus one lone-pair donor (furan, thiophene, pyrrole):
///   O or S with only its two ring bonds, or N with only single bonds; the four carbons
///   each have one double bond (Hückel 6 pi).
/// A ring atom's double bond may lie in an already aromatic fused ring (naphthalene,
/// indole). Each ring lists its atoms in ascending order; rings are sorted, so the result
/// does not depend on how the SMILES was written.
[[nodiscard]] std::vector<std::vector<std::size_t>> aromatic_rings(const Molecule& molecule);

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

/// Which atoms give the same NMR signal. Topological symmetry (symmetry_classes) is not
/// enough: geometry, slow rotation and stereochemistry make some topologically equivalent
/// groups inequivalent in NMR. Because ARISE must never make an "expert" call that could
/// rule out the true structure, uncertain cases give a range between two partitions.
/// Classes are numbered from 0 in order of first appearance, as in symmetry_classes.
struct NmrEquivalence {
    /// Splits that always happen at room temperature (the fewest signals):
    /// - E/Z: two identical groups on one end of a non-aromatic C=C or C=N whose other end
    ///   carries two different groups (a nitrogen lone pair counts as a group), e.g. the
    ///   =C(CH3)2 methyls of linalool, the alpha-CH2 groups of an oxime;
    /// - slow rotation: the two identical N-groups of an amide or thioamide (DMF).
    std::vector<std::size_t> fewest;
    /// Additionally, splits that may or may not show (the most signals):
    /// - borderline rotation: identical N-groups of carbamates, ureas and enaminones;
    /// - when any stereocentre exists (sp3 C with four different groups, sulfoxide S,
    ///   P with three different groups) or at least two prochiral centres:
    ///   diastereotopic pairs (an sp3 C or Si, or a quaternary N+, carrying two identical
    ///   groups and two other different ones, e.g. valine's isopropyl; a P, S, Se or I
    ///   centre with one identical pair among three or more groups, e.g. ArI(OAc)2), and
    ///   every stereocentre on its own, since
    ///   unknown relative configuration can break a symmetry (mannitol has 3 carbon
    ///   environments, sorbitol 6). Over-splitting only widens the range.
    std::vector<std::size_t> most;
};

[[nodiscard]] NmrEquivalence nmr_equivalence(const Molecule& molecule);

/// Number of distinct classes among the carbon atoms of `classes` (one entry per atom).
[[nodiscard]] int carbon_environments(const Molecule& molecule,
                                      const std::vector<std::size_t>& classes);

/// Functional groups a chemist can confirm or rule out: fine-grained motifs that NMR can
/// tell apart, as signed off in docs/chem_group_definitions.md (each maps to a ChEBI
/// class). Counts are per carbonyl carbon, per OH, per N, per bond or per ring.
enum class Group : std::uint8_t {
    CarboxylicAcid,       // C(=O)OH, acyl side C or H
    Aldehyde,             // C(=O)H, carbonyl C with no single bond to O, N or S
    Ketone,               // C(=O) single-bonded to exactly two carbons
    CarboxylicEster,      // C(=O)-O-C, acyl side C or H; lactones, formates
    CarbonateEster,       // C-O-C(=O)-O-C
    CarboxylicAnhydride,  // C(=O)-O-C(=O), counted per bridging O
    Carboxamide,          // C(=O)-N, acyl side C or H; lactams, formamides
    CarbamateEster,       // C-O-C(=O)-N, or C=N-O-C(=O)-N (O-carbamoyl oximes)
    Urea,                 // N-C(=O)-N
    Alcohol,              // OH on an sp3 carbon
    Enol,                 // OH on a non-aromatic C=C carbon
    Phenol,               // OH on a carbon of an all-carbon aromatic ring
    Pyridone,             // 2-/4-pyridone ring, or OH on a pyridine carbon
    Ether,                // C-O-C, neither carbon a carbonyl C, O not in an aromatic ring
    Amine,                // N with single bonds to C/H only, not to C=O/C=S, not aromatic
    Sulfonamide,          // S(=O)2-N, counted per S
    Hydrazine,            // N-N single bond outside aromatic rings
    Hydroxylamine,        // N-O single bond, N with only single bonds
    Imine,                // C=N, N bonded to no O or N
    Oxime,                // C=N-OH
    Nitrile,              // C#N
    Alkene,               // C=C outside aromatic rings
    BenzeneRing,          // aromatic six-membered all-carbon ring
    PyridineRing,         // aromatic six-membered ring with one N
    FuranRing,            // aromatic five-membered ring with one O
    PyrroleRing,          // aromatic five-membered ring with one N
    ThiopheneRing,        // aromatic five-membered ring with one S
};

/// All groups, in enum order.
inline constexpr std::array kAllGroups{
    Group::CarboxylicAcid,  Group::Aldehyde,       Group::Ketone,
    Group::CarboxylicEster, Group::CarbonateEster, Group::CarboxylicAnhydride,
    Group::Carboxamide,     Group::CarbamateEster, Group::Urea,
    Group::Alcohol,         Group::Enol,           Group::Phenol,
    Group::Pyridone,        Group::Ether,          Group::Amine,
    Group::Sulfonamide,     Group::Hydrazine,      Group::Hydroxylamine,
    Group::Imine,           Group::Oxime,          Group::Nitrile,
    Group::Alkene,          Group::BenzeneRing,    Group::PyridineRing,
    Group::FuranRing,       Group::PyrroleRing,    Group::ThiopheneRing,
};

/// "carboxylic acid", "aldehyde", ...
[[nodiscard]] std::string_view to_string(Group group) noexcept;

/// The atoms carrying each NMR hallmark of `group`, as named in the "Key atom" column of
/// docs/chem_group_definitions.md: "C=O carbon"; an ester also "O-C carbon"; "C-OH
/// carbon", "ring carbons", ... Empty for groups whose hallmark is not on carbon
/// (sulfonamide, hydrazine, hydroxylamine).
[[nodiscard]] std::span<const std::string_view> hallmarks(Group group) noexcept;

/// One occurrence of a group (one unit of count_group) and, per hallmark of
/// hallmarks(group), the atoms carrying it, in ascending order (e.g. a ketone: {{the C=O
/// carbon}}; an ester: {{C=O carbon}, {O-C carbon}}; a benzene ring: {{its six carbons}}).
struct GroupOccurrence {
    std::vector<std::vector<std::size_t>> key_atoms;
};

/// Every occurrence of `group`, with its key atoms.
[[nodiscard]] std::vector<GroupOccurrence> find_group(const Molecule& molecule, Group group);

/// How many times `group` occurs in the molecule (per carbonyl carbon, per OH, per C=C,
/// per ring, ...): find_group(molecule, group).size().
[[nodiscard]] int count_group(const Molecule& molecule, Group group);

}  // namespace arise::chem
