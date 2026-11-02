#include "arise/chem/perception.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <set>
#include <utility>

namespace arise::chem {
namespace {

constexpr int kAromaticBondCode = 4;     // distinct from bond orders 1-3
constexpr int kDelocalizedBondCode = 5;  // a bond to a resonance-shared terminal atom
constexpr int kDelocalizedCharge = 100;  // invariant of a resonance-shared terminal atom

// A ring's atoms in walking order.
using Ring = std::vector<std::size_t>;

bool is_ring_element(Element element) noexcept {
    return element == Element::C || element == Element::N || element == Element::O ||
           element == Element::S;
}

// Index into molecule.bonds() of the bond between a and b, if any.
std::optional<std::size_t> find_bond(const Molecule& molecule, std::size_t a,
                                     std::size_t b) noexcept {
    for (const auto bond : molecule.bonds_of(a)) {
        if (molecule.neighbor(bond, a) == b) {
            return bond;
        }
    }
    return std::nullopt;
}

BondOrder order_of(const Molecule& molecule, std::size_t bond) noexcept {
    return molecule.bonds()[bond].order;
}

// Every simple cycle of `size` atoms among C/N/O/S, each listed once in walking order,
// starting at its smallest atom.
class RingFinder {
public:
    RingFinder(const Molecule& molecule, std::size_t size) : molecule_{&molecule}, size_{size} {}

    std::vector<Ring> find() {
        for (std::size_t start = 0; start < molecule_->atoms().size(); ++start) {
            if (is_ring_element(molecule_->atoms()[start].element)) {
                path_ = {start};
                extend(start);
            }
        }
        return std::move(rings_);
    }

private:
    void extend(std::size_t start) {
        const std::size_t last = path_.back();
        for (const auto bond : molecule_->bonds_of(last)) {
            const std::size_t next = molecule_->neighbor(bond, last);
            if (path_.size() == size_) {
                if (next == start) {
                    record();
                }
                continue;
            }
            if (next <= start || !is_ring_element(molecule_->atoms()[next].element) ||
                std::find(path_.begin(), path_.end(), next) != path_.end()) {
                continue;
            }
            path_.push_back(next);
            extend(start);
            path_.pop_back();
        }
    }

    // Each cycle is walked in both directions; keep the one whose second atom is smaller.
    void record() {
        if (path_[1] < path_.back()) {
            rings_.push_back(path_);
        }
    }

    const Molecule* molecule_;  // not null
    std::size_t size_;
    std::vector<std::size_t> path_;
    std::vector<Ring> rings_;
};

// The single double-bond partner of `atom`, or nothing when the atom has no double bond,
// several, or a triple bond.
std::optional<std::size_t> sole_double_partner(const Molecule& molecule, std::size_t atom) {
    std::optional<std::size_t> partner;
    for (const auto bond : molecule.bonds_of(atom)) {
        const auto order = order_of(molecule, bond);
        if (order == BondOrder::Triple || (order == BondOrder::Double && partner.has_value())) {
            return std::nullopt;
        }
        if (order == BondOrder::Double) {
            partner = molecule.neighbor(bond, atom);
        }
    }
    return partner;
}

bool has_only_single_bonds(const Molecule& molecule, std::size_t atom) {
    const auto bonds = molecule.bonds_of(atom);
    return std::all_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return order_of(molecule, bond) == BondOrder::Single;
    });
}

// Whether ring atom `i` has exactly one double bond, lying in the ring or in an aromatic
// ring already found (fused systems).
bool pi_carbon(const Molecule& molecule, const Ring& ring, std::size_t i,
               const std::set<std::size_t>& aromatic_bonds) {
    const std::size_t atom = ring[i];
    const auto partner = sole_double_partner(molecule, atom);
    if (!partner.has_value()) {
        return false;
    }
    const std::size_t n = ring.size();
    const bool in_ring = *partner == ring[(i + 1) % n] || *partner == ring[(i + n - 1) % n];
    const auto bond = find_bond(molecule, atom, *partner);
    return in_ring || (bond.has_value() && aromatic_bonds.contains(*bond));
}

// Six-membered: every atom C or N with one pi bond. Five-membered: one lone-pair donor
// (O or S with only its two ring bonds, or N with only single bonds) and four carbons with
// one pi bond each (Huckel 6 pi: furan, thiophene, pyrrole).
bool is_aromatic(const Molecule& molecule, const Ring& ring,
                 const std::set<std::size_t>& aromatic_bonds) {
    const auto& atoms = molecule.atoms();
    if (ring.size() == 6) {
        for (std::size_t i = 0; i < 6; ++i) {
            const auto element = atoms[ring[i]].element;
            if ((element != Element::C && element != Element::N) ||
                !pi_carbon(molecule, ring, i, aromatic_bonds)) {
                return false;
            }
        }
        return true;
    }
    int donors = 0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto atom = ring[i];
        const auto element = atoms[atom].element;
        if (element == Element::C) {
            if (!pi_carbon(molecule, ring, i, aromatic_bonds)) {
                return false;
            }
            continue;
        }
        const bool donor = has_only_single_bonds(molecule, atom) &&
                           (element == Element::N || molecule.bonds_of(atom).size() == 2);
        if (!donor || (element != Element::N && element != Element::O && element != Element::S)) {
            return false;
        }
        ++donors;
    }
    return donors == 1;
}

struct Aromaticity {
    std::vector<Ring> rings;      // in walking order
    std::set<std::size_t> bonds;  // bond indices inside aromatic rings
    std::vector<bool> atoms;      // per atom: in an aromatic ring
};

Aromaticity perceive_aromaticity(const Molecule& molecule) {
    Aromaticity result;
    result.atoms.assign(molecule.atoms().size(), false);
    auto pending = RingFinder{molecule, 6}.find();
    auto five = RingFinder{molecule, 5}.find();
    pending.insert(pending.end(), five.begin(), five.end());
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto it = pending.begin(); it != pending.end();) {
            if (!is_aromatic(molecule, *it, result.bonds)) {
                ++it;
                continue;
            }
            const std::size_t n = it->size();
            for (std::size_t i = 0; i < n; ++i) {
                if (const auto bond = find_bond(molecule, (*it)[i], (*it)[(i + 1) % n])) {
                    result.bonds.insert(*bond);
                }
                result.atoms[(*it)[i]] = true;
            }
            result.rings.push_back(*it);
            it = pending.erase(it);
            changed = true;
        }
    }
    return result;
}

// Dense ranks of `signatures`: equal signatures get equal ranks, ordered by value, so
// the ranks do not depend on atom order.
std::vector<int> rank(const std::vector<std::vector<int>>& signatures) {
    std::map<std::vector<int>, int> ranks;
    for (const auto& signature : signatures) {
        ranks.emplace(signature, 0);
    }
    int next = 0;
    for (auto& [signature, value] : ranks) {
        value = next++;
    }
    std::vector<int> result;
    result.reserve(signatures.size());
    for (const auto& signature : signatures) {
        result.push_back(ranks.at(signature));
    }
    return result;
}

std::size_t distinct(const std::vector<int>& labels) {
    return std::set<int>(labels.begin(), labels.end()).size();
}

bool is_carbonyl_carbon(const Molecule& molecule, std::size_t atom) {
    // A neutral C=O: an oxonium C=[O+] (an ylide form PubChem sometimes lists) is not one.
    if (molecule.atoms()[atom].element != Element::C || molecule.atoms()[atom].charge != 0) {
        return false;
    }
    const auto bonds = molecule.bonds_of(atom);
    return std::any_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        const auto& other = molecule.atoms()[molecule.neighbor(bond, atom)];
        return order_of(molecule, bond) == BondOrder::Double && other.element == Element::O &&
               other.charge == 0;
    });
}

// ---- Symmetry refinement ----------------------------------------------------------------

bool is_heavy(const Molecule& molecule, std::size_t atom) {
    return molecule.atoms()[atom].element != Element::H;
}

// Bonds from a centre atom to terminal atoms that share a charge by resonance: terminal
// heavy neighbours of one element, without hydrogens, each either double-bonded and
// neutral or single-bonded with charge -1, at least one of each kind (the O of nitro,
// carboxylate and sulfonate groups: [N+](=O)[O-], C(=O)[O-], S(=O)(=O)[O-]).
std::set<std::size_t> delocalized_bonds(const Molecule& molecule) {
    const auto& atoms = molecule.atoms();
    std::set<std::size_t> result;
    for (std::size_t centre = 0; centre < atoms.size(); ++centre) {
        std::map<Element, std::vector<std::size_t>> by_element;
        std::map<Element, std::pair<bool, bool>> kinds;  // {double neutral, single anion}
        for (const auto bond : molecule.bonds_of(centre)) {
            const auto end = molecule.neighbor(bond, centre);
            const auto& atom = atoms[end];
            if (atom.element == Element::H || molecule.hydrogen_count(end) != 0 ||
                molecule.bonds_of(end).size() != 1) {
                continue;
            }
            const auto order = order_of(molecule, bond);
            const bool double_neutral = order == BondOrder::Double && atom.charge == 0;
            const bool single_anion = order == BondOrder::Single && atom.charge == -1;
            if (double_neutral || single_anion) {
                by_element[atom.element].push_back(bond);
                kinds[atom.element].first = kinds[atom.element].first || double_neutral;
                kinds[atom.element].second = kinds[atom.element].second || single_anion;
            }
        }
        for (const auto& [element, bonds] : by_element) {
            if (kinds[element].first && kinds[element].second) {
                result.insert(bonds.begin(), bonds.end());
            }
        }
    }
    return result;
}

// Starting invariants: element, isotope, charge, total hydrogens, heavy neighbours. Atoms
// that share a charge by resonance (delocalized_bonds) get one charge invariant.
// Hydrogen atoms enter only through hydrogen_count, never as neighbours, so a methyl
// written with explicit [H] atoms matches one with implicit hydrogens.
std::vector<int> initial_labels(const Molecule& molecule) {
    const auto& atoms = molecule.atoms();
    std::vector<bool> shared(atoms.size(), false);
    for (const auto bond : delocalized_bonds(molecule)) {
        shared[molecule.bonds()[bond].first] = true;
        shared[molecule.bonds()[bond].second] = true;
    }
    std::vector<std::vector<int>> signatures(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto bonds = molecule.bonds_of(i);
        const auto heavy_degree = std::count_if(bonds.begin(), bonds.end(), [&](std::size_t b) {
            return is_heavy(molecule, molecule.neighbor(b, i));
        });
        signatures[i] = {
            static_cast<int>(atoms[i].element),
            atoms[i].isotope,
            shared[i] && molecule.bonds_of(i).size() == 1 ? kDelocalizedCharge : atoms[i].charge,
            molecule.hydrogen_count(i),
            static_cast<int>(heavy_degree),
        };
    }
    return rank(signatures);
}

// Refines `labels` on neighbour labels until the number of classes stops growing. Bonds
// of aromatic rings count as one kind, so Kekule benzene stays symmetric; so do the bonds
// to resonance-shared terminal atoms, so the two O of a nitro group stay equivalent.
std::vector<int> refine(const Molecule& molecule, const Aromaticity& aromaticity,
                        std::vector<int> labels) {
    const std::size_t n = molecule.atoms().size();
    const auto delocalized = delocalized_bonds(molecule);
    std::vector<std::vector<int>> signatures(n);
    for (std::size_t round = 0; round <= n; ++round) {
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<std::pair<int, int>> around;
            for (const auto bond : molecule.bonds_of(i)) {
                const auto other = molecule.neighbor(bond, i);
                if (!is_heavy(molecule, other)) {
                    continue;
                }
                int code = static_cast<int>(order_of(molecule, bond));
                if (aromaticity.bonds.contains(bond)) {
                    code = kAromaticBondCode;
                } else if (delocalized.contains(bond)) {
                    code = kDelocalizedBondCode;
                }
                around.emplace_back(code, labels[other]);
            }
            std::sort(around.begin(), around.end());
            signatures[i] = {labels[i]};
            for (const auto& [code, label] : around) {
                signatures[i].push_back(code);
                signatures[i].push_back(label);
            }
        }
        auto refined = rank(signatures);
        const bool stable = distinct(refined) == distinct(labels);
        labels = std::move(refined);
        if (stable) {
            break;
        }
    }
    return labels;
}

std::vector<std::size_t> first_appearance(const std::vector<int>& labels) {
    std::map<int, std::size_t> first_seen;
    std::vector<std::size_t> classes;
    classes.reserve(labels.size());
    for (const int label : labels) {
        const auto [it, inserted] = first_seen.emplace(label, first_seen.size());
        classes.push_back(it->second);
    }
    return classes;
}

// ---- NMR equivalence rules ----------------------------------------------------------------

// An atom to separate from its class. Picks of one class normally share their new label,
// so symmetric copies stay paired; `own_label` gives the atom a label of its own.
struct Pick {
    std::size_t atom;
    bool own_label = false;
};

// Heavy neighbours of `atom`, except `skip`.
std::vector<std::size_t> heavy_neighbors(const Molecule& molecule, std::size_t atom,
                                         std::optional<std::size_t> skip = std::nullopt) {
    std::vector<std::size_t> result;
    for (const auto bond : molecule.bonds_of(atom)) {
        const auto other = molecule.neighbor(bond, atom);
        if (is_heavy(molecule, other) && other != skip) {
            result.push_back(other);
        }
    }
    return result;
}

bool has_double_bond_to(const Molecule& molecule, std::size_t atom, Element element) {
    const auto bonds = molecule.bonds_of(atom);
    return std::any_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return order_of(molecule, bond) == BondOrder::Double &&
               molecule.atoms()[molecule.neighbor(bond, atom)].element == element;
    });
}

// The lower-index member of a pair of identical groups, if `groups` is exactly such a pair.
std::optional<std::size_t> identical_pair(const std::vector<std::size_t>& groups,
                                          const std::vector<int>& labels) {
    if (groups.size() != 2 || labels[groups[0]] != labels[groups[1]]) {
        return std::nullopt;
    }
    return std::min(groups[0], groups[1]);
}

// Whether the substituents of `atom` other than `skip` are two *different* groups. A
// nitrogen's lone pair counts as a group, so a C=N nitrogen always qualifies.
bool two_different_groups(const Molecule& molecule, std::size_t atom, std::size_t skip,
                          const std::vector<int>& labels) {
    const auto heavy = heavy_neighbors(molecule, atom, skip);
    const int hydrogens = molecule.hydrogen_count(atom);
    if (molecule.atoms()[atom].element == Element::N) {
        return heavy.size() + static_cast<std::size_t>(hydrogens) == 1;
    }
    if (heavy.size() + static_cast<std::size_t>(hydrogens) != 2) {
        return false;
    }
    if (heavy.size() == 2) {
        return labels[heavy[0]] != labels[heavy[1]];
    }
    return heavy.size() == 1;  // one group and one H; two H are identical
}

// E/Z: identical groups on one end of a non-aromatic C=C or C=N, different on the other.
std::vector<Pick> ez_splits(const Molecule& molecule, const Aromaticity& aromaticity,
                            const std::vector<int>& labels) {
    std::vector<Pick> picks;
    const auto& atoms = molecule.atoms();
    for (std::size_t b = 0; b < molecule.bonds().size(); ++b) {
        const auto& bond = molecule.bonds()[b];
        if (bond.order != BondOrder::Double || aromaticity.bonds.contains(b)) {
            continue;
        }
        for (const auto& [end, other] :
             {std::pair{bond.first, bond.second}, std::pair{bond.second, bond.first}}) {
            const auto other_element = atoms[other].element;
            if (atoms[end].element != Element::C ||
                (other_element != Element::C && other_element != Element::N)) {
                continue;
            }
            const auto pair = identical_pair(heavy_neighbors(molecule, end, other), labels);
            if (pair.has_value() && molecule.hydrogen_count(end) == 0 &&
                two_different_groups(molecule, other, end, labels)) {
                picks.push_back(Pick{.atom = *pair});
            }
        }
    }
    return picks;
}

enum class Barrier : std::uint8_t { Slow, Borderline };

// The C=O/C=S (or vinylogous) partner that slows rotation about N's bond to it, with how
// slow: amides and thioamides are slow; carbamates, ureas and enaminones borderline.
std::optional<std::pair<std::size_t, Barrier>> rotation_partner(const Molecule& molecule,
                                                                std::size_t nitrogen) {
    const auto& atoms = molecule.atoms();
    for (const auto k : heavy_neighbors(molecule, nitrogen)) {
        if (atoms[k].element != Element::C) {
            continue;
        }
        if (has_double_bond_to(molecule, k, Element::O) ||
            has_double_bond_to(molecule, k, Element::S)) {
            const auto others = heavy_neighbors(molecule, k, nitrogen);
            const bool hetero = std::any_of(others.begin(), others.end(), [&](std::size_t o) {
                const auto e = atoms[o].element;
                return (e == Element::N || e == Element::O) &&
                       find_bond(molecule, k, o).has_value() &&
                       order_of(molecule, *find_bond(molecule, k, o)) == BondOrder::Single;
            });
            return std::pair{k, hetero ? Barrier::Borderline : Barrier::Slow};
        }
        // Enaminone: N-C=C-C(=O).
        for (const auto bond : molecule.bonds_of(k)) {
            const auto c2 = molecule.neighbor(bond, k);
            if (order_of(molecule, bond) != BondOrder::Double || atoms[c2].element != Element::C) {
                continue;
            }
            for (const auto c3 : heavy_neighbors(molecule, c2, k)) {
                if (atoms[c3].element == Element::C &&
                    has_double_bond_to(molecule, c3, Element::O)) {
                    return std::pair{k, Barrier::Borderline};
                }
            }
        }
    }
    return std::nullopt;
}

// Identical groups on an amide-type nitrogen, for the given barrier class.
std::vector<Pick> rotation_splits(const Molecule& molecule, const std::vector<int>& labels,
                                  Barrier barrier) {
    std::vector<Pick> picks;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (molecule.atoms()[i].element != Element::N || !has_only_single_bonds(molecule, i) ||
            molecule.hydrogen_count(i) != 0) {
            continue;
        }
        const auto partner = rotation_partner(molecule, i);
        if (!partner.has_value() || partner->second != barrier) {
            continue;
        }
        if (const auto pair =
                identical_pair(heavy_neighbors(molecule, i, partner->first), labels)) {
            picks.push_back(Pick{.atom = *pair});
        }
    }
    return picks;
}

// Whether all substituent groups of `atom` (heavy neighbours by label, plus each H) differ.
bool all_groups_differ(const Molecule& molecule, std::size_t atom, const std::vector<int>& labels) {
    const auto heavy = heavy_neighbors(molecule, atom);
    std::set<int> seen;
    for (const auto n : heavy) {
        if (!seen.insert(labels[n]).second) {
            return false;
        }
    }
    return molecule.hydrogen_count(atom) <= 1;
}

bool is_stereocentre(const Molecule& molecule, std::size_t atom, const std::vector<int>& labels) {
    const auto element = molecule.atoms()[atom].element;
    const auto groups = heavy_neighbors(molecule, atom).size() +
                        static_cast<std::size_t>(molecule.hydrogen_count(atom));
    if (element == Element::C) {
        return has_only_single_bonds(molecule, atom) && groups == 4 &&
               all_groups_differ(molecule, atom, labels);
    }
    if (element == Element::S) {  // sulfoxide: S(=O) with two different groups
        return groups == 3 && has_double_bond_to(molecule, atom, Element::O) &&
               molecule.bonds_of(atom).size() == 3 && all_groups_differ(molecule, atom, labels);
    }
    if (element == Element::P) {  // three different groups, with or without P=O/P=S
        return groups >= 3 && all_groups_differ(molecule, atom, labels);
    }
    return false;
}

// Whether `atom` can be a prochiral centre: an sp3 C or Si, or a quaternary N+, with four
// groups (Caabc); or a P, S, Se or I centre with at least three groups (ArI(OAc)2, the two
// O of a sulfone). Neutral amines are left out: nitrogen inversion averages their groups.
bool can_be_prochiral(const Molecule& molecule, std::size_t atom, std::size_t groups) {
    const auto& a = molecule.atoms()[atom];
    switch (a.element) {
        case Element::C:
        case Element::Si:
            return has_only_single_bonds(molecule, atom) && groups == 4;
        case Element::N:
            return a.charge == 1 && has_only_single_bonds(molecule, atom) && groups == 4;
        case Element::P:
        case Element::S:
        case Element::Se:
        case Element::I:
            return groups >= 3;
        default:
            return false;
    }
}

// The lower-index member of the identical pair on a prochiral centre (can_be_prochiral):
// two identical groups and the others all different, at most one H.
std::optional<std::size_t> prochiral_pair(const Molecule& molecule, std::size_t atom,
                                          const std::vector<int>& labels) {
    const auto heavy = heavy_neighbors(molecule, atom);
    const int hydrogens = molecule.hydrogen_count(atom);
    const auto groups = heavy.size() + static_cast<std::size_t>(hydrogens);
    if (!can_be_prochiral(molecule, atom, groups) || hydrogens > 1) {
        return std::nullopt;
    }
    std::map<int, std::vector<std::size_t>> by_label;
    for (const auto h : heavy) {
        by_label[labels[h]].push_back(h);
    }
    std::optional<std::size_t> pair;
    for (const auto& [label, members] : by_label) {
        if (members.size() > 2 || (members.size() == 2 && pair.has_value())) {
            return std::nullopt;  // three identical groups, or two identical pairs
        }
        if (members.size() == 2) {
            pair = std::min(members[0], members[1]);
        }
    }
    const bool rest_differ = by_label.size() + static_cast<std::size_t>(hydrogens) == groups - 1;
    return rest_differ ? pair : std::nullopt;
}

// Diastereotopic splits, which may or may not show: when the molecule has a stereocentre
// or at least two prochiral centres (these make each other's groups diastereotopic), one
// of each prochiral pair, and every stereocentre on its own (unknown configuration can
// break a symmetry). Over-splitting only widens the range, the safe direction.
std::vector<Pick> stereo_splits(const Molecule& molecule, const std::vector<int>& labels) {
    std::vector<Pick> stereocentres;
    std::vector<Pick> prochiral;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (is_stereocentre(molecule, i, labels)) {
            stereocentres.push_back(Pick{.atom = i, .own_label = true});
        } else if (const auto pair = prochiral_pair(molecule, i, labels)) {
            prochiral.push_back(Pick{.atom = *pair});
        }
    }
    if (stereocentres.empty() && prochiral.size() < 2) {
        return {};
    }
    stereocentres.insert(stereocentres.end(), prochiral.begin(), prochiral.end());
    return stereocentres;
}

// Gives each class of picked atoms one new label and refines, repeating until the rules
// pick nothing new. Picked atoms of one class share a label, so symmetric copies (both
// "cis" methyls of a symmetric molecule) stay equivalent to each other.
template <typename Rules>
std::vector<int> split_until_stable(const Molecule& molecule, const Aromaticity& aromaticity,
                                    std::vector<int> labels, Rules rules) {
    for (std::size_t round = 0; round <= molecule.atoms().size(); ++round) {
        const auto picks = rules(labels);
        if (picks.empty()) {
            break;
        }
        const std::size_t before = distinct(labels);
        int next = *std::max_element(labels.begin(), labels.end()) + 1;
        std::map<int, int> new_label;
        for (const auto& pick : picks) {
            if (pick.own_label) {
                labels[pick.atom] = next++;
                continue;
            }
            const auto [it, inserted] = new_label.emplace(labels[pick.atom], next);
            next += inserted ? 1 : 0;
            labels[pick.atom] = it->second;
        }
        labels = refine(molecule, aromaticity, std::move(labels));
        if (distinct(labels) == before) {
            break;
        }
    }
    return labels;
}

// ---- Functional groups (docs/chem_group_definitions.md) -----------------------------------

bool is_element(const Molecule& molecule, std::size_t atom, Element element) {
    return molecule.atoms()[atom].element == element;
}

// What a carbonyl carbon's single-bonded neighbours are, to tell its family apart.
struct CarbonylSides {
    std::vector<std::size_t> oxygens;    // single-bonded O
    std::vector<std::size_t> nitrogens;  // single-bonded N
    std::vector<std::size_t> carbons;    // single-bonded C
    bool other = false;                  // single-bonded S, halogen, ...
};

CarbonylSides carbonyl_sides(const Molecule& molecule, std::size_t carbon) {
    CarbonylSides sides;
    for (const auto bond : molecule.bonds_of(carbon)) {
        const auto n = molecule.neighbor(bond, carbon);
        if (order_of(molecule, bond) != BondOrder::Single || !is_heavy(molecule, n)) {
            continue;
        }
        if (is_element(molecule, n, Element::O)) {
            sides.oxygens.push_back(n);
        } else if (is_element(molecule, n, Element::N)) {
            sides.nitrogens.push_back(n);
        } else if (is_element(molecule, n, Element::C)) {
            sides.carbons.push_back(n);
        } else {
            sides.other = true;
        }
    }
    return sides;
}

// An OH, or an O- with a single heavy neighbour, which counts as the OH of its neutral
// parent (carboxylate as acid, alkoxide as alcohol; user decision 2026-11-01).
bool is_hydroxyl_like(const Molecule& molecule, std::size_t oxygen) {
    if (molecule.hydrogen_count(oxygen) > 0) {
        return true;
    }
    return molecule.atoms()[oxygen].charge == -1 && heavy_neighbors(molecule, oxygen).size() == 1;
}

// An O bonded by single bonds to `carbonyl` and to another carbon that is not a carbonyl C.
bool ester_oxygen(const Molecule& molecule, std::size_t oxygen, std::size_t carbonyl) {
    if (molecule.hydrogen_count(oxygen) > 0) {
        return false;
    }
    const auto others = heavy_neighbors(molecule, oxygen, carbonyl);
    return others.size() == 1 && is_element(molecule, others[0], Element::C) &&
           !is_carbonyl_carbon(molecule, others[0]);
}

// A carbamate's O bonded to an oxime N (C=N-O-C(=O)-N: methomyl, aldicarb).
bool oxime_carbamate_oxygen(const Molecule& molecule, std::size_t oxygen, std::size_t carbonyl) {
    const auto others = heavy_neighbors(molecule, oxygen, carbonyl);
    return molecule.hydrogen_count(oxygen) == 0 && others.size() == 1 &&
           is_element(molecule, others[0], Element::N) &&
           has_double_bond_to(molecule, others[0], Element::C);
}

bool carbonyl_is(const Molecule& molecule, std::size_t carbon, Group group) {
    const auto sides = carbonyl_sides(molecule, carbon);
    const std::size_t o = sides.oxygens.size();
    const std::size_t n = sides.nitrogens.size();
    if (sides.other) {
        return false;
    }
    const auto any_ester_o = [&] {
        return std::any_of(sides.oxygens.begin(), sides.oxygens.end(),
                           [&](std::size_t x) { return ester_oxygen(molecule, x, carbon); });
    };
    switch (group) {
        case Group::CarboxylicAcid:
            return o == 1 && n == 0 && is_hydroxyl_like(molecule, sides.oxygens[0]);
        case Group::Aldehyde:
            return o == 0 && n == 0 && molecule.hydrogen_count(carbon) > 0;
        case Group::Ketone:
            return o == 0 && n == 0 && sides.carbons.size() == 2;
        case Group::CarboxylicEster:
            return o == 1 && n == 0 && any_ester_o();
        case Group::CarbonateEster:
            return o == 2 && n == 0 && ester_oxygen(molecule, sides.oxygens[0], carbon) &&
                   ester_oxygen(molecule, sides.oxygens[1], carbon);
        case Group::Carboxamide:
            return o == 0 && n == 1;
        case Group::CarbamateEster:
            return o == 1 && n == 1 &&
                   (any_ester_o() || oxime_carbamate_oxygen(molecule, sides.oxygens[0], carbon));
        case Group::Urea:
            return o == 0 && n == 2;
        default:
            return false;
    }
}

// The single heavy neighbour of a hydroxyl O, if the O is an OH with exactly one. By the
// neutral parent (user decision 2026-11-01), an O- with one heavy neighbour counts as the
// OH it came from: carboxylate, alkoxide, phenoxide.
std::optional<std::size_t> hydroxyl_carbon(const Molecule& molecule, std::size_t oxygen) {
    const auto heavy = heavy_neighbors(molecule, oxygen);
    if (!is_hydroxyl_like(molecule, oxygen) || heavy.size() != 1 ||
        !is_element(molecule, heavy[0], Element::C)) {
        return std::nullopt;
    }
    return heavy[0];
}

bool has_nonaromatic_double_bond(const Molecule& molecule, const Aromaticity& aromaticity,
                                 std::size_t atom, Element partner) {
    const auto bonds = molecule.bonds_of(atom);
    return std::any_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return order_of(molecule, bond) == BondOrder::Double && !aromaticity.bonds.contains(bond) &&
               is_element(molecule, molecule.neighbor(bond, atom), partner);
    });
}

// Carbon atoms in an aromatic ring made only of carbons (benzene rings).
std::vector<bool> carbocyclic_aromatic_atoms(const Molecule& molecule,
                                             const Aromaticity& aromaticity) {
    std::vector<bool> result(molecule.atoms().size(), false);
    for (const auto& ring : aromaticity.rings) {
        const bool all_carbon = std::all_of(ring.begin(), ring.end(), [&](std::size_t a) {
            return is_element(molecule, a, Element::C);
        });
        for (const auto a : ring) {
            result[a] = result[a] || all_carbon;
        }
    }
    return result;
}

// Whether an aromatic ring has `size` atoms, exactly one of them `hetero` (or none when
// `hetero` is C), the rest carbon.
bool ring_is(const Molecule& molecule, const Ring& ring, std::size_t size, Element hetero) {
    if (ring.size() != size) {
        return false;
    }
    const auto count = std::count_if(
        ring.begin(), ring.end(), [&](std::size_t a) { return is_element(molecule, a, hetero); });
    const auto carbons = std::count_if(ring.begin(), ring.end(), [&](std::size_t a) {
        return is_element(molecule, a, Element::C);
    });
    if (hetero == Element::C) {
        return std::cmp_equal(carbons, size);
    }
    return count == 1 && std::cmp_equal(carbons + 1, size);
}

using Occurrences = std::vector<GroupOccurrence>;

// One occurrence whose hallmarks are carried by `atoms` (one list per hallmark), sorted.
GroupOccurrence occurrence(std::vector<std::vector<std::size_t>> atoms) {
    for (auto& list : atoms) {
        std::sort(list.begin(), list.end());
    }
    return GroupOccurrence{.key_atoms = std::move(atoms)};
}

// The ester O-C carbons of an ester carbonyl carbon.
std::vector<std::size_t> ester_carbons(const Molecule& molecule, std::size_t carbonyl) {
    std::vector<std::size_t> carbons;
    for (const auto o : carbonyl_sides(molecule, carbonyl).oxygens) {
        if (ester_oxygen(molecule, o, carbonyl)) {
            const auto others = heavy_neighbors(molecule, o, carbonyl);
            carbons.insert(carbons.end(), others.begin(), others.end());
        }
    }
    return carbons;
}

// Carbonyl-carbon groups (key atom the C=O carbon; an ester also its O-C carbon) and
// anhydrides (one per bridging O; key atoms both C=O carbons).
Occurrences find_carbonyl_group(const Molecule& molecule, Group group) {
    Occurrences found;
    const std::size_t n = molecule.atoms().size();
    if (group == Group::CarboxylicAnhydride) {
        for (std::size_t o = 0; o < n; ++o) {
            const auto heavy = heavy_neighbors(molecule, o);
            if (is_element(molecule, o, Element::O) && heavy.size() == 2 &&
                is_carbonyl_carbon(molecule, heavy[0]) && is_carbonyl_carbon(molecule, heavy[1])) {
                found.push_back(occurrence({heavy}));
            }
        }
        return found;
    }
    for (std::size_t c = 0; c < n; ++c) {
        if (!is_carbonyl_carbon(molecule, c) || !carbonyl_is(molecule, c, group)) {
            continue;
        }
        if (group == Group::CarboxylicEster) {
            found.push_back(occurrence({{c}, ester_carbons(molecule, c)}));
        } else {
            found.push_back(occurrence({{c}}));
        }
    }
    return found;
}

Occurrences find_oxygen_group(const Molecule& molecule, const Aromaticity& aromaticity,
                              Group group) {
    const auto carbocyclic = carbocyclic_aromatic_atoms(molecule, aromaticity);
    Occurrences found;
    for (std::size_t o = 0; o < molecule.atoms().size(); ++o) {
        if (!is_element(molecule, o, Element::O) || aromaticity.atoms[o]) {
            continue;
        }
        const auto carbon = hydroxyl_carbon(molecule, o);
        bool match = false;
        switch (group) {
            case Group::Alcohol:
                match = carbon.has_value() && has_only_single_bonds(molecule, *carbon);
                break;
            case Group::Enol:
                match = carbon.has_value() &&
                        has_nonaromatic_double_bond(molecule, aromaticity, *carbon, Element::C);
                break;
            case Group::Phenol:
                match = carbon.has_value() && carbocyclic[*carbon];
                break;
            case Group::Ether: {
                const auto heavy = heavy_neighbors(molecule, o);
                if (molecule.hydrogen_count(o) == 0 && heavy.size() == 2 &&
                    has_only_single_bonds(molecule, o) &&
                    std::all_of(heavy.begin(), heavy.end(), [&](std::size_t c) {
                        return is_element(molecule, c, Element::C) &&
                               !is_carbonyl_carbon(molecule, c);
                    })) {
                    found.push_back(occurrence({heavy}));
                }
                break;
            }
            default:
                break;
        }
        if (match) {
            found.push_back(occurrence({{*carbon}}));
        }
    }
    return found;
}

// A 2-/4-pyridone ring: five C and one N with only single bonds, a C=O on a ring carbon,
// and two C=C in the ring.
bool is_pyridone_ring(const Molecule& molecule, const Ring& ring) {
    int nitrogens = 0;
    int carbonyls = 0;
    int ring_double = 0;
    for (std::size_t i = 0; i < 6; ++i) {
        const auto a = ring[i];
        nitrogens +=
            is_element(molecule, a, Element::N) && has_only_single_bonds(molecule, a) ? 1 : 0;
        carbonyls += is_carbonyl_carbon(molecule, a) ? 1 : 0;
        const auto bond = find_bond(molecule, a, ring[(i + 1) % 6]);
        ring_double += bond.has_value() && order_of(molecule, *bond) == BondOrder::Double ? 1 : 0;
    }
    const auto carbons = std::count_if(ring.begin(), ring.end(), [&](std::size_t a) {
        return is_element(molecule, a, Element::C);
    });
    return nitrogens == 1 && carbonyls == 1 && ring_double == 2 && carbons == 5;
}

// Pyridone rings (key atom: the ring C=O carbon) plus OH groups on aromatic pyridine
// carbons (the hydroxypyridine tautomer; key atom: the C-OH carbon).
Occurrences find_pyridones(const Molecule& molecule, const Aromaticity& aromaticity) {
    Occurrences found;
    for (const auto& ring : RingFinder{molecule, 6}.find()) {
        if (is_pyridone_ring(molecule, ring)) {
            const auto carbonyl = std::find_if(ring.begin(), ring.end(), [&](std::size_t a) {
                return is_carbonyl_carbon(molecule, a);
            });
            found.push_back(occurrence({{*carbonyl}}));
        }
    }
    for (std::size_t o = 0; o < molecule.atoms().size(); ++o) {
        if (!is_element(molecule, o, Element::O)) {
            continue;
        }
        const auto carbon = hydroxyl_carbon(molecule, o);
        if (!carbon.has_value() || !aromaticity.atoms[*carbon]) {
            continue;
        }
        const bool in_pyridine =
            std::any_of(aromaticity.rings.begin(), aromaticity.rings.end(), [&](const Ring& ring) {
                return ring_is(molecule, ring, 6, Element::N) &&
                       std::find(ring.begin(), ring.end(), *carbon) != ring.end();
            });
        if (in_pyridine) {
            found.push_back(occurrence({{*carbon}}));
        }
    }
    return found;
}

bool amine_nitrogen(const Molecule& molecule, const Aromaticity& aromaticity, std::size_t n) {
    if (!is_element(molecule, n, Element::N) || aromaticity.atoms[n] ||
        !has_only_single_bonds(molecule, n)) {
        return false;
    }
    const auto heavy = heavy_neighbors(molecule, n);
    return std::all_of(heavy.begin(), heavy.end(), [&](std::size_t c) {
        return is_element(molecule, c, Element::C) && !is_carbonyl_carbon(molecule, c) &&
               !has_double_bond_to(molecule, c, Element::S);
    });
}

// Sulfonamide S atoms; the hallmark is not on carbon, so no key atoms.
Occurrences find_sulfonamides(const Molecule& molecule) {
    Occurrences found;
    for (std::size_t s = 0; s < molecule.atoms().size(); ++s) {
        if (!is_element(molecule, s, Element::S)) {
            continue;
        }
        int oxo = 0;
        bool nitrogen = false;
        for (const auto bond : molecule.bonds_of(s)) {
            const auto other = molecule.neighbor(bond, s);
            oxo += order_of(molecule, bond) == BondOrder::Double &&
                           is_element(molecule, other, Element::O)
                       ? 1
                       : 0;
            nitrogen = nitrogen || (order_of(molecule, bond) == BondOrder::Single &&
                                    is_element(molecule, other, Element::N));
        }
        if (oxo == 2 && nitrogen) {
            found.push_back(occurrence({}));
        }
    }
    return found;
}

// Bonds between `a` and `b` atoms of `order`, outside aromatic rings, passing `accept`,
// as (the `a` atom, the `b` atom).
template <typename Accept>
std::vector<std::pair<std::size_t, std::size_t>> find_bonds(const Molecule& molecule,
                                                            const Aromaticity& aromaticity,
                                                            BondOrder order, Element a, Element b,
                                                            Accept accept) {
    std::vector<std::pair<std::size_t, std::size_t>> found;
    for (std::size_t i = 0; i < molecule.bonds().size(); ++i) {
        const auto& bond = molecule.bonds()[i];
        if (bond.order != order || aromaticity.bonds.contains(i)) {
            continue;
        }
        const auto ea = molecule.atoms()[bond.first].element;
        const auto eb = molecule.atoms()[bond.second].element;
        if (ea == a && eb == b && accept(bond.first, bond.second)) {
            found.emplace_back(bond.first, bond.second);
        } else if (ea == b && eb == a && accept(bond.second, bond.first)) {
            found.emplace_back(bond.second, bond.first);
        }
    }
    return found;
}

// For a C=N nitrogen: whether it carries a single-bonded O (oxime) or N (hydrazone).
bool nitrogen_has_single(const Molecule& molecule, std::size_t n, Element element,
                         bool with_hydrogen) {
    const auto bonds = molecule.bonds_of(n);
    return std::any_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        const auto other = molecule.neighbor(bond, n);
        return order_of(molecule, bond) == BondOrder::Single &&
               is_element(molecule, other, element) &&
               (!with_hydrogen || molecule.hydrogen_count(other) > 0);
    });
}

std::vector<std::pair<std::size_t, std::size_t>> nitrogen_bond_group_bonds(
    const Molecule& molecule, const Aromaticity& aromaticity, Group group) {
    const auto always = [](std::size_t, std::size_t) { return true; };
    switch (group) {
        case Group::Hydrazine:
            return find_bonds(molecule, aromaticity, BondOrder::Single, Element::N, Element::N,
                              [&](std::size_t a, std::size_t b) {
                                  return !aromaticity.atoms[a] && !aromaticity.atoms[b];
                              });
        case Group::Hydroxylamine:
            return find_bonds(molecule, aromaticity, BondOrder::Single, Element::N, Element::O,
                              [&](std::size_t n, std::size_t o) {
                                  // An N-oxide's N+-O- is not a hydroxylamine (user
                                  // decision 2026-11-01).
                                  const bool n_oxide = molecule.atoms()[n].charge == 1 &&
                                                       molecule.atoms()[o].charge == -1;
                                  return has_only_single_bonds(molecule, n) && !n_oxide &&
                                         !aromaticity.atoms[n] && !aromaticity.atoms[o];
                              });
        case Group::Imine:
            return find_bonds(molecule, aromaticity, BondOrder::Double, Element::C, Element::N,
                              [&](std::size_t, std::size_t n) {
                                  return molecule.atoms()[n].charge == 0 &&
                                         !nitrogen_has_single(molecule, n, Element::O, false) &&
                                         !nitrogen_has_single(molecule, n, Element::N, false);
                              });
        case Group::Oxime:
            return find_bonds(molecule, aromaticity, BondOrder::Double, Element::C, Element::N,
                              [&](std::size_t, std::size_t n) {
                                  return molecule.atoms()[n].charge == 0 &&
                                         nitrogen_has_single(molecule, n, Element::O, true);
                              });
        case Group::Nitrile:  // not an isocyanide, [C-]#[N+]
            return find_bonds(molecule, aromaticity, BondOrder::Triple, Element::C, Element::N,
                              [&](std::size_t c, std::size_t n) {
                                  return molecule.atoms()[c].charge == 0 &&
                                         molecule.atoms()[n].charge == 0;
                              });
        case Group::Alkene:
            return find_bonds(molecule, aromaticity, BondOrder::Double, Element::C, Element::C,
                              always);
        default:
            return {};
    }
}

// The bond groups: key atom the C of C=N and C#N, both C of C=C, none for N-N and N-O.
Occurrences find_nitrogen_bond_group(const Molecule& molecule, const Aromaticity& aromaticity,
                                     Group group) {
    Occurrences found;
    for (const auto& [a, b] : nitrogen_bond_group_bonds(molecule, aromaticity, group)) {
        if (group == Group::Hydrazine || group == Group::Hydroxylamine) {
            found.push_back(occurrence({}));
        } else if (group == Group::Alkene) {
            found.push_back(occurrence({{a, b}}));
        } else {
            found.push_back(occurrence({{a}}));  // the carbon of C=N or C#N
        }
    }
    return found;
}

// Aromatic rings of `size` with `hetero` (C for all-carbon); key atoms the ring carbons.
Occurrences find_rings(const Molecule& molecule, const Aromaticity& aromaticity, std::size_t size,
                       Element hetero) {
    Occurrences found;
    for (const auto& ring : aromaticity.rings) {
        if (ring_is(molecule, ring, size, hetero)) {
            std::vector<std::size_t> carbons;
            std::copy_if(ring.begin(), ring.end(), std::back_inserter(carbons),
                         [&](std::size_t a) { return is_element(molecule, a, Element::C); });
            found.push_back(occurrence({carbons}));
        }
    }
    return found;
}

}  // namespace

std::optional<nmr1d::CarbonType> carbon_type(const Molecule& molecule, std::size_t atom) noexcept {
    if (molecule.atoms()[atom].element != Element::C) {
        return std::nullopt;
    }
    switch (molecule.hydrogen_count(atom)) {
        case 0:
            return nmr1d::CarbonType::C;
        case 1:
            return nmr1d::CarbonType::CH;
        case 2:
            return nmr1d::CarbonType::CH2;
        case 3:
            return nmr1d::CarbonType::CH3;
        default:
            return std::nullopt;
    }
}

std::vector<std::vector<std::size_t>> aromatic_rings(const Molecule& molecule) {
    auto rings = perceive_aromaticity(molecule).rings;
    for (auto& ring : rings) {
        std::sort(ring.begin(), ring.end());
    }
    std::sort(rings.begin(), rings.end());
    return rings;
}

std::vector<std::size_t> symmetry_classes(const Molecule& molecule) {
    const auto aromaticity = perceive_aromaticity(molecule);
    return first_appearance(refine(molecule, aromaticity, initial_labels(molecule)));
}

NmrEquivalence nmr_equivalence(const Molecule& molecule) {
    const auto aromaticity = perceive_aromaticity(molecule);
    const auto topological = refine(molecule, aromaticity, initial_labels(molecule));
    const auto hard = [&](const std::vector<int>& labels) {
        auto picks = ez_splits(molecule, aromaticity, labels);
        auto rotation = rotation_splits(molecule, labels, Barrier::Slow);
        picks.insert(picks.end(), rotation.begin(), rotation.end());
        return picks;
    };
    const auto all = [&](const std::vector<int>& labels) {
        auto picks = hard(labels);
        auto borderline = rotation_splits(molecule, labels, Barrier::Borderline);
        picks.insert(picks.end(), borderline.begin(), borderline.end());
        auto stereo = stereo_splits(molecule, labels);
        picks.insert(picks.end(), stereo.begin(), stereo.end());
        return picks;
    };
    const auto fewest = split_until_stable(molecule, aromaticity, topological, hard);
    const auto most = split_until_stable(molecule, aromaticity, fewest, all);
    return NmrEquivalence{.fewest = first_appearance(fewest), .most = first_appearance(most)};
}

int carbon_environments(const Molecule& molecule, const std::vector<std::size_t>& classes) {
    std::set<std::size_t> carbon_classes;
    for (std::size_t i = 0; i < molecule.atoms().size(); ++i) {
        if (molecule.atoms()[i].element == Element::C) {
            carbon_classes.insert(classes[i]);
        }
    }
    return static_cast<int>(carbon_classes.size());
}

std::string_view to_string(Group group) noexcept {
    switch (group) {
        case Group::CarboxylicAcid:
            return "carboxylic acid";
        case Group::Aldehyde:
            return "aldehyde";
        case Group::Ketone:
            return "ketone";
        case Group::CarboxylicEster:
            return "carboxylic ester";
        case Group::CarbonateEster:
            return "carbonate ester";
        case Group::CarboxylicAnhydride:
            return "carboxylic anhydride";
        case Group::Carboxamide:
            return "carboxamide";
        case Group::CarbamateEster:
            return "carbamate ester";
        case Group::Urea:
            return "urea";
        case Group::Alcohol:
            return "alcohol";
        case Group::Enol:
            return "enol";
        case Group::Phenol:
            return "phenol";
        case Group::Pyridone:
            return "pyridone";
        case Group::Ether:
            return "ether";
        case Group::Amine:
            return "amine";
        case Group::Sulfonamide:
            return "sulfonamide";
        case Group::Hydrazine:
            return "hydrazine";
        case Group::Hydroxylamine:
            return "hydroxylamine";
        case Group::Imine:
            return "imine";
        case Group::Oxime:
            return "oxime";
        case Group::Nitrile:
            return "nitrile";
        case Group::Alkene:
            return "alkene";
        case Group::BenzeneRing:
            return "benzene ring";
        case Group::PyridineRing:
            return "pyridine ring";
        case Group::FuranRing:
            return "furan ring";
        case Group::PyrroleRing:
            return "pyrrole ring";
        case Group::ThiopheneRing:
            return "thiophene ring";
    }
    return "?";
}

std::span<const std::string_view> hallmarks(Group group) noexcept {
    static constexpr std::array<std::string_view, 1> kCarbonyl{"C=O carbon"};
    static constexpr std::array<std::string_view, 2> kEster{"C=O carbon", "O-C carbon"};
    static constexpr std::array<std::string_view, 1> kHydroxyl{"C-OH carbon"};
    static constexpr std::array<std::string_view, 1> kPyridone{"C=O or C-OH carbon"};
    static constexpr std::array<std::string_view, 1> kEther{"O-C carbon"};
    static constexpr std::array<std::string_view, 1> kAmine{"N-C carbon"};
    static constexpr std::array<std::string_view, 1> kImine{"C=N carbon"};
    static constexpr std::array<std::string_view, 1> kNitrile{"C#N carbon"};
    static constexpr std::array<std::string_view, 1> kAlkene{"C=C carbon"};
    static constexpr std::array<std::string_view, 1> kRing{"ring carbons"};
    switch (group) {
        case Group::CarboxylicAcid:
        case Group::Aldehyde:
        case Group::Ketone:
        case Group::CarbonateEster:
        case Group::CarboxylicAnhydride:
        case Group::Carboxamide:
        case Group::CarbamateEster:
        case Group::Urea:
            return kCarbonyl;
        case Group::CarboxylicEster:
            return kEster;
        case Group::Alcohol:
        case Group::Enol:
        case Group::Phenol:
            return kHydroxyl;
        case Group::Pyridone:
            return kPyridone;
        case Group::Ether:
            return kEther;
        case Group::Amine:
            return kAmine;
        case Group::Imine:
        case Group::Oxime:
            return kImine;
        case Group::Nitrile:
            return kNitrile;
        case Group::Alkene:
            return kAlkene;
        case Group::BenzeneRing:
        case Group::PyridineRing:
        case Group::FuranRing:
        case Group::PyrroleRing:
        case Group::ThiopheneRing:
            return kRing;
        case Group::Sulfonamide:
        case Group::Hydrazine:
        case Group::Hydroxylamine:
            return {};
    }
    return {};
}

std::vector<GroupOccurrence> find_group(const Molecule& molecule, Group group) {
    const auto aromaticity = perceive_aromaticity(molecule);
    switch (group) {
        case Group::CarboxylicAcid:
        case Group::Aldehyde:
        case Group::Ketone:
        case Group::CarboxylicEster:
        case Group::CarbonateEster:
        case Group::CarboxylicAnhydride:
        case Group::Carboxamide:
        case Group::CarbamateEster:
        case Group::Urea:
            return find_carbonyl_group(molecule, group);
        case Group::Alcohol:
        case Group::Enol:
        case Group::Phenol:
        case Group::Ether:
            return find_oxygen_group(molecule, aromaticity, group);
        case Group::Pyridone:
            return find_pyridones(molecule, aromaticity);
        case Group::Amine: {
            Occurrences found;
            for (std::size_t a = 0; a < molecule.atoms().size(); ++a) {
                if (amine_nitrogen(molecule, aromaticity, a)) {
                    found.push_back(occurrence({heavy_neighbors(molecule, a)}));
                }
            }
            return found;
        }
        case Group::Sulfonamide:
            return find_sulfonamides(molecule);
        case Group::Hydrazine:
        case Group::Hydroxylamine:
        case Group::Imine:
        case Group::Oxime:
        case Group::Nitrile:
        case Group::Alkene:
            return find_nitrogen_bond_group(molecule, aromaticity, group);
        case Group::BenzeneRing:
            return find_rings(molecule, aromaticity, 6, Element::C);
        case Group::PyridineRing:
            return find_rings(molecule, aromaticity, 6, Element::N);
        case Group::FuranRing:
            return find_rings(molecule, aromaticity, 5, Element::O);
        case Group::PyrroleRing:
            return find_rings(molecule, aromaticity, 5, Element::N);
        case Group::ThiopheneRing:
            return find_rings(molecule, aromaticity, 5, Element::S);
    }
    return {};
}

int count_group(const Molecule& molecule, Group group) {
    return static_cast<int>(find_group(molecule, group).size());
}

}  // namespace arise::chem
