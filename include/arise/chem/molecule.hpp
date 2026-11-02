#pragma once

#include "arise/core/result.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace arise::chem {

/// Elements the reader accepts: those found in organic NMR samples (the same set as
/// core::Formula).
enum class Element : std::uint8_t { H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I };

/// "C", "Cl", ...
[[nodiscard]] std::string_view symbol(Element element) noexcept;

struct Atom {
    Element element = Element::C;
    std::int8_t charge = 0;
    /// Hydrogens carried by this atom: implicit for unbracketed atoms, as written for
    /// bracket atoms. Explicit [H] atoms are separate atoms (see Molecule::hydrogen_count).
    std::uint8_t hydrogens = 0;
    /// Mass number, 0 for natural abundance.
    std::uint16_t isotope = 0;
    /// Written as a bracket atom, e.g. [CH] or [13CH4].
    bool bracket = false;

    friend bool operator==(const Atom&, const Atom&) = default;
};

enum class BondOrder : std::uint8_t { Single = 1, Double = 2, Triple = 3 };

struct Bond {
    std::size_t first = 0;
    std::size_t second = 0;
    BondOrder order = BondOrder::Single;

    friend bool operator==(const Bond&, const Bond&) = default;
};

/// One connected molecule as an atom/bond graph. Atoms are numbered in SMILES order.
class Molecule {
public:
    [[nodiscard]] std::span<const Atom> atoms() const noexcept { return atoms_; }
    [[nodiscard]] std::span<const Bond> bonds() const noexcept { return bonds_; }

    /// Indices into bonds() of the bonds touching `atom`, in SMILES order.
    [[nodiscard]] std::span<const std::size_t> bonds_of(std::size_t atom) const noexcept;

    /// The atom at the other end of bond `bond` from `atom`.
    [[nodiscard]] std::size_t neighbor(std::size_t bond, std::size_t atom) const noexcept;

    /// Hydrogens on `atom`: its own count plus explicit [H] atoms bonded to it.
    [[nodiscard]] int hydrogen_count(std::size_t atom) const noexcept;

    /// Bond orders to all neighbours plus the atom's own hydrogens.
    [[nodiscard]] int valence(std::size_t atom) const noexcept;

private:
    friend core::Result<Molecule> parse_smiles(std::string_view smiles);

    std::vector<Atom> atoms_;
    std::vector<Bond> bonds_;
    // Compressed adjacency: the bonds of atom i are adjacency_[offsets_[i] .. offsets_[i+1]).
    std::vector<std::size_t> offsets_;
    std::vector<std::size_t> adjacency_;
};

/// Reads a SMILES string in Kekulé form, as PubChem writes it, into one molecule:
/// - Unbracketed organic-subset atoms B C N O P S F Cl Br I get implicit hydrogens up to
///   their lowest normal valence that fits (OpenSMILES): B 3; C 4; N 3, 5; O 2; P 3, 5;
///   S 2, 4, 6; F, Cl, Br 1; I 1, 3, 5 (beyond OpenSMILES: hypervalent iodine).
/// - Bracket atoms: isotope, element (also H, Si, Se), hydrogen count, charge;
///   stereo marks (@, @@) are accepted and ignored.
/// - Bonds - = # (and directional / \ as single), branches, ring closures 0-9 and %nn
///   with the bond order on either end.
/// ParseError (naming the 1-based column) for empty input, lowercase aromatic atoms,
/// unsupported elements or notation, more than one component ('.'), unbalanced
/// branches, unclosed or conflicting ring closures, duplicate bonds and dangling bonds.
[[nodiscard]] core::Result<Molecule> parse_smiles(std::string_view smiles);

/// Why a parsed molecule is not a usable candidate.
enum class StructureIssue : std::uint8_t {
    Charged,         // an ion: formal charges that do not sum to zero
    Radical,         // unpaired electrons: valence below or between normal valences
    UnusualValence,  // valence above the atom's highest normal valence
};

/// "charged", "radical" or "unusual valence".
[[nodiscard]] std::string_view to_string(StructureIssue issue) noexcept;

/// The first issue found, checked in this order: net charge, radical, unusual valence.
/// Charge-separated neutral molecules are fine: nitro groups, N-oxides, zwitterions and
/// isocyanides, which PubChem writes with formal charges ([N+](=O)[O-]).
/// Normal valences: H, F, Cl, Br 1; I 1, 3, 5; O 2; B 3; C, Si 4; N, P 3, 5; S, Se 2, 4, 6.
/// A charge shifts them: groups 15-17 gain a bond when positive and lose one when negative
/// (N+ 4, O- 1, S+ 3/5/7), boron the reverse (B- 4), carbon and silicon make 3 either way;
/// second-period atoms shift only their lowest valence. Neutral N at valence 5 with only
/// single bonds is an unusual valence (five sigma bonds).
/// Empty for a closed-shell net-neutral molecule.
[[nodiscard]] std::optional<StructureIssue> structure_issue(const Molecule& molecule);

}  // namespace arise::chem
