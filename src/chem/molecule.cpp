#include "arise/chem/molecule.hpp"

#include "arise/core/error.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <string>
#include <utility>

namespace arise::chem {
namespace {

struct ElementInfo {
    Element element;
    std::string_view symbol;
    int atomic_number;
    int group;                    // periodic table group (1-18)
    std::array<int, 3> valences;  // normal valences, ascending; unused slots are 0
    bool organic_subset;          // may be written without brackets
};

// Indexed by Element. Normal valences follow OpenSMILES (organic subset); Si and Se, which
// OpenSMILES leaves to bracket atoms, are checked against RDKit by the valence probes
// (tests/data/reference/valence_probes.tsv).
constexpr std::array kElements{
    ElementInfo{
        .element = Element::H,
        .symbol = "H",
        .atomic_number = 1,
        .group = 1,
        .valences = {1, 0, 0},
        .organic_subset = false,
    },
    ElementInfo{
        .element = Element::B,
        .symbol = "B",
        .atomic_number = 5,
        .group = 13,
        .valences = {3, 0, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::C,
        .symbol = "C",
        .atomic_number = 6,
        .group = 14,
        .valences = {4, 0, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::N,
        .symbol = "N",
        .atomic_number = 7,
        .group = 15,
        .valences = {3, 5, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::O,
        .symbol = "O",
        .atomic_number = 8,
        .group = 16,
        .valences = {2, 0, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::F,
        .symbol = "F",
        .atomic_number = 9,
        .group = 17,
        .valences = {1, 0, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::Si,
        .symbol = "Si",
        .atomic_number = 14,
        .group = 14,
        .valences = {4, 0, 0},
        .organic_subset = false,
    },
    ElementInfo{
        .element = Element::P,
        .symbol = "P",
        .atomic_number = 15,
        .group = 15,
        .valences = {3, 5, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::S,
        .symbol = "S",
        .atomic_number = 16,
        .group = 16,
        .valences = {2, 4, 6},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::Cl,
        .symbol = "Cl",
        .atomic_number = 17,
        .group = 17,
        .valences = {1, 0, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::Se,
        .symbol = "Se",
        .atomic_number = 34,
        .group = 16,
        .valences = {2, 4, 6},
        .organic_subset = false,
    },
    ElementInfo{
        .element = Element::Br,
        .symbol = "Br",
        .atomic_number = 35,
        .group = 17,
        .valences = {1, 0, 0},
        .organic_subset = true,
    },
    ElementInfo{
        .element = Element::I,
        .symbol = "I",
        .atomic_number = 53,
        .group = 17,
        .valences = {1, 3, 5},  // hypervalent iodine reagents (user decision 2026-11-01)
        .organic_subset = true,
    },
};
constexpr std::size_t kMaxRings = 100;
constexpr int kLastSecondPeriod = 10;  // neon
constexpr int kBoronGroup = 13;
constexpr int kCarbonGroup = 14;
constexpr int kNitrogenGroup = 15;

const ElementInfo& info(Element element) noexcept {
    return kElements[static_cast<std::size_t>(element)];
}

std::optional<Element> find_element(std::string_view symbol) noexcept {
    for (const auto& element : kElements) {
        if (element.symbol == symbol) {
            return element.element;
        }
    }
    return std::nullopt;
}

using Valences = std::array<int, 3>;

// Normal valences of an atom. A charge changes how many bonds an atom makes: in groups
// 15-17 a positive charge adds one and a negative one removes one (N+ 4, O- 1, S- 1/3/5,
// I+ 2/4/6); boron gains one when negative and loses one when positive (B- 4, B+ 2);
// carbon and silicon make three either way (C+ 3, C- 3). Second-period atoms get only
// their lowest valence shifted (N+ 4, not 6). The valence probes check all of this
// against RDKit (tests/data/reference/valence_probes.tsv).
Valences normal_valences(const Atom& atom) noexcept {
    const auto& element = info(atom.element);
    if (atom.charge == 0) {
        return element.valences;
    }
    Valences valences = element.atomic_number <= kLastSecondPeriod
                            ? Valences{element.valences[0], 0, 0}
                            : element.valences;
    const int charge = +atom.charge;  // unary plus promotes the int8_t to int
    int shift = 0;
    if (element.group == kBoronGroup) {
        shift = -charge;
    } else if (element.group == kCarbonGroup) {
        shift = -(charge < 0 ? -charge : charge);
    } else if (element.group >= kNitrogenGroup) {
        shift = charge;
    } else {
        return {0, 0, 0};  // charged hydrogen: not a closed-shell organic atom
    }
    for (auto& valence : valences) {
        valence = valence != 0 && valence + shift > 0 ? valence + shift : 0;
    }
    return valences;
}

int max_valence(const Valences& valences) noexcept {
    return *std::max_element(valences.begin(), valences.end());
}

bool is_normal_valence(const Valences& valences, int valence) noexcept {
    return std::any_of(valences.begin(), valences.end(),
                       [valence](int normal) { return normal != 0 && normal == valence; });
}

bool is_upper(char c) noexcept {
    return c >= 'A' && c <= 'Z';
}
bool is_lower(char c) noexcept {
    return c >= 'a' && c <= 'z';
}
bool is_digit(char c) noexcept {
    return c >= '0' && c <= '9';
}

core::Error error_at(std::size_t column, const std::string& what) {
    return core::Error{
        .code = core::ErrorCode::ParseError,
        .message = "column " + std::to_string(column) + ": " + what,
    };
}

std::string quoted(char c) {
    return "'" + std::string(1, c) + "'";
}

std::optional<BondOrder> bond_symbol(char c) noexcept {
    switch (c) {
        case '-':
        case '/':  // directional single bonds: stereo is ignored
        case '\\':
            return BondOrder::Single;
        case '=':
            return BondOrder::Double;
        case '#':
            return BondOrder::Triple;
        default:
            return std::nullopt;
    }
}

struct Graph {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
};

struct PendingBond {
    BondOrder order;
    std::size_t column;
};

// The bond order written at one end of a ring closure, if any. A plain flag rather than
// std::optional: GCC 13 at -O3 reports a false "maybe-uninitialized" for the optional.
struct RingBondEnd {
    BondOrder order = BondOrder::Single;
    bool written = false;
};

struct OpenRing {
    std::size_t atom;
    RingBondEnd end;
};

struct OpenBranch {
    std::size_t atom;
    std::size_t atoms_before;
    std::size_t column;
};

// Reads one SMILES string, token by token.
class Parser {
public:
    explicit Parser(std::string_view text) : text_{text} {}

    core::Result<Graph> run() {
        if (text_.empty()) {
            return core::Error{.code = core::ErrorCode::ParseError, .message = "empty SMILES"};
        }
        while (pos_ < text_.size()) {
            auto step = next();
            if (!step.ok()) {
                return std::move(step).error();
            }
        }
        auto finished = finish();
        if (!finished.ok()) {
            return std::move(finished).error();
        }
        return std::move(graph_);
    }

private:
    [[nodiscard]] std::size_t column() const noexcept { return pos_ + 1; }

    core::Result<void> next() {
        const char c = text_[pos_];
        if (is_upper(c)) {
            return organic_atom();
        }
        if (c == '[') {
            return bracket_atom();
        }
        if (const auto order = bond_symbol(c)) {
            return bond(*order);
        }
        if (is_digit(c) || c == '%') {
            return ring_closure();
        }
        if (c == '(') {
            return open_branch();
        }
        if (c == ')') {
            return close_branch();
        }
        return unexpected(c);
    }

    core::Result<void> unexpected(char c) const {
        if (is_lower(c)) {
            return error_at(column(),
                            "aromatic (lowercase) atoms are not supported yet; write "
                            "the SMILES in Kekulé form (explicit double bonds)");
        }
        if (c == ':') {
            return error_at(column(),
                            "aromatic bonds (':') are not supported yet; write the "
                            "SMILES in Kekulé form");
        }
        if (c == '.') {
            return error_at(column(), "more than one component ('.'); expected a single molecule");
        }
        if (c == '$') {
            return error_at(column(), "'$' (quadruple bond) is not supported");
        }
        return error_at(column(), "unexpected character " + quoted(c));
    }

    core::Result<void> organic_atom() {
        const std::size_t start = column();
        const auto two = text_.substr(pos_, 2);
        if (two == "Cl" || two == "Br") {
            const auto halogen = find_element(two);
            assert(halogen.has_value());  // both are in kElements
            pos_ += 2;
            return add_atom(Atom{.element = *halogen});
        }
        if (two.size() == 2 && is_lower(two[1]) && find_element(two).has_value()) {
            return error_at(start, "'" + std::string{two} + "' must be written in brackets ([" +
                                       std::string{two} + "])");
        }
        const auto element = find_element(text_.substr(pos_, 1));
        if (!element.has_value()) {
            return error_at(start, "unsupported element " + quoted(text_[pos_]));
        }
        if (!info(*element).organic_subset) {
            return error_at(
                start, "'" + std::string{info(*element).symbol} + "' must be written in brackets");
        }
        ++pos_;
        return add_atom(Atom{.element = *element});
    }

    core::Result<void> bracket_atom() {
        const std::size_t open = column();
        const auto close = text_.find(']', pos_);
        if (close == std::string_view::npos) {
            return error_at(open, "bracket atom is missing its ']'");
        }
        if (close == pos_ + 1) {
            return error_at(open, "empty bracket atom");
        }
        BracketReader reader{text_, pos_ + 1, close};
        auto atom = reader.read();
        if (!atom.ok()) {
            return std::move(atom).error();
        }
        pos_ = close + 1;
        return add_atom(atom.value());
    }

    // The inside of one bracket atom: [isotope] symbol [@|@@] [H[n]] [charge].
    class BracketReader {
    public:
        BracketReader(std::string_view text, std::size_t begin, std::size_t end)
            : text_{text}, pos_{begin}, end_{end} {}

        core::Result<Atom> read() {
            Atom atom{.bracket = true};
            atom.isotope = static_cast<std::uint16_t>(number(999));
            auto element = symbol();
            if (!element.ok()) {
                return std::move(element).error();
            }
            atom.element = element.value();
            auto stereo = skip_stereo();
            if (!stereo.ok()) {
                return std::move(stereo).error();
            }
            if (pos_ < end_ && text_[pos_] == 'H') {
                ++pos_;
                const bool has_count = pos_ < end_ && is_digit(text_[pos_]);
                atom.hydrogens = static_cast<std::uint8_t>(has_count ? number(9) : 1);
            }
            atom.charge = static_cast<std::int8_t>(charge());
            if (pos_ < end_ && text_[pos_] == ':') {
                return error_at(pos_ + 1, "atom classes (':n') are not supported");
            }
            if (pos_ < end_) {
                return error_at(pos_ + 1, "unexpected " + quoted(text_[pos_]) + " in bracket atom");
            }
            return atom;
        }

    private:
        // Digits at the cursor, capped at `limit`; 0 when there are none.
        int number(int limit) {
            int value = 0;
            while (pos_ < end_ && is_digit(text_[pos_])) {
                value = (value * 10) + (text_[pos_] - '0');
                value = value > limit ? limit : value;
                ++pos_;
            }
            return value;
        }

        core::Result<Element> symbol() {
            const std::size_t start = pos_ + 1;
            if (pos_ >= end_ || !is_upper(text_[pos_])) {
                if (pos_ < end_ && is_lower(text_[pos_])) {
                    return error_at(start,
                                    "aromatic (lowercase) atoms are not supported yet; "
                                    "write the SMILES in Kekulé form");
                }
                return error_at(start, "bracket atom needs an element symbol");
            }
            const std::size_t length = (pos_ + 1 < end_ && is_lower(text_[pos_ + 1])) ? 2 : 1;
            const auto text = text_.substr(pos_, length);
            const auto element = find_element(text);
            if (!element.has_value()) {
                return error_at(start, "unsupported element '" + std::string{text} + "'");
            }
            pos_ += length;
            return *element;
        }

        core::Result<void> skip_stereo() {
            const std::size_t start = pos_ + 1;
            int marks = 0;
            while (pos_ < end_ && text_[pos_] == '@') {
                ++marks;
                ++pos_;
            }
            if (marks > 2 ||
                (marks > 0 && pos_ < end_ && is_upper(text_[pos_]) && text_[pos_] != 'H')) {
                return error_at(start, "only @ and @@ stereo marks are supported");
            }
            return {};
        }

        int charge() {
            if (pos_ >= end_ || (text_[pos_] != '+' && text_[pos_] != '-')) {
                return 0;
            }
            const char sign = text_[pos_];
            const int unit = sign == '+' ? 1 : -1;
            ++pos_;
            if (pos_ < end_ && is_digit(text_[pos_])) {
                return unit * number(9);
            }
            int count = 1;
            while (pos_ < end_ && text_[pos_] == sign) {
                ++count;
                ++pos_;
            }
            return unit * count;
        }

        std::string_view text_;
        std::size_t pos_;
        std::size_t end_;
    };

    core::Result<void> add_atom(const Atom& atom) {
        graph_.atoms.push_back(atom);
        const std::size_t index = graph_.atoms.size() - 1;
        if (previous_.has_value()) {
            const auto order = pending_.has_value() ? pending_->order : BondOrder::Single;
            auto added = add_bond(*previous_, index, order, column());
            if (!added.ok()) {
                return added;
            }
        }
        pending_.reset();
        previous_ = index;
        return {};
    }

    core::Result<void> bond(BondOrder order) {
        if (!previous_.has_value()) {
            return error_at(column(), "a bond must follow an atom");
        }
        if (pending_.has_value()) {
            return error_at(column(), "two bond symbols in a row");
        }
        pending_ = PendingBond{.order = order, .column = column()};
        ++pos_;
        return {};
    }

    core::Result<void> ring_closure() {
        const std::size_t start = column();
        std::size_t ring = 0;
        if (text_[pos_] == '%') {
            if (pos_ + 2 >= text_.size() || !is_digit(text_[pos_ + 1]) ||
                !is_digit(text_[pos_ + 2])) {
                return error_at(start, "'%' must be followed by two digits");
            }
            const int tens = text_[pos_ + 1] - '0';
            const int units = text_[pos_ + 2] - '0';
            const int two_digits = (tens * 10) + units;
            ring = static_cast<std::size_t>(two_digits);
            pos_ += 3;
        } else {
            const int digit = text_[pos_] - '0';
            ring = static_cast<std::size_t>(digit);
            ++pos_;
        }
        if (!previous_.has_value()) {
            return error_at(start, "a ring closure must follow an atom");
        }
        RingBondEnd here;
        if (pending_.has_value()) {
            here = RingBondEnd{.order = pending_->order, .written = true};
        }
        pending_.reset();
        auto& open = rings_[ring];
        if (!open.has_value()) {
            open = OpenRing{.atom = *previous_, .end = here};
            ring_columns_[ring] = start;
            return {};
        }
        const auto opened = *open;
        open.reset();
        if (opened.atom == *previous_) {
            return error_at(start, "ring " + std::to_string(ring) + " closes on itself");
        }
        if (opened.end.written && here.written && opened.end.order != here.order) {
            return error_at(start, "ring " + std::to_string(ring) + " bond orders conflict");
        }
        const auto order = opened.end.written ? opened.end.order : here.order;
        return add_bond(opened.atom, *previous_, order, start);
    }

    core::Result<void> open_branch() {
        if (!previous_.has_value()) {
            return error_at(column(), "a branch must follow an atom");
        }
        if (pending_.has_value()) {
            return error_at(pending_->column, "a bond must be inside the branch, e.g. C(=O)");
        }
        branches_.push_back(OpenBranch{
            .atom = *previous_,
            .atoms_before = graph_.atoms.size(),
            .column = column(),
        });
        ++pos_;
        return {};
    }

    core::Result<void> close_branch() {
        if (branches_.empty()) {
            return error_at(column(), "')' closes no branch");
        }
        if (pending_.has_value()) {
            return error_at(pending_->column, "the bond is not followed by an atom");
        }
        if (graph_.atoms.size() == branches_.back().atoms_before) {
            return error_at(column(), "empty branch");
        }
        previous_ = branches_.back().atom;
        branches_.pop_back();
        ++pos_;
        return {};
    }

    // Bonds are stored with first < second, so one comparison finds duplicates. Chain
    // bonds already come in that order; a ring closed after a branch may not ("CC(C1C)1").
    core::Result<void> add_bond(std::size_t a, std::size_t b, BondOrder order, std::size_t at) {
        if (a > b) {
            std::swap(a, b);
        }
        for (const auto& existing : graph_.bonds) {
            if (existing.first == a && existing.second == b) {
                return error_at(at, "duplicate bond between the same two atoms");
            }
        }
        graph_.bonds.push_back(Bond{.first = a, .second = b, .order = order});
        return {};
    }

    core::Result<void> finish() const {
        if (pending_.has_value()) {
            return error_at(pending_->column, "the bond is not followed by an atom");
        }
        if (!branches_.empty()) {
            return error_at(branches_.back().column, "branch '(' is not closed");
        }
        for (std::size_t ring = 0; ring < kMaxRings; ++ring) {
            if (rings_[ring].has_value()) {
                return error_at(ring_columns_[ring],
                                "ring " + std::to_string(ring) + " is not closed");
            }
        }
        return {};
    }

    std::string_view text_;
    std::size_t pos_ = 0;
    Graph graph_;
    std::optional<std::size_t> previous_;
    std::optional<PendingBond> pending_;
    std::vector<OpenBranch> branches_;
    std::array<std::optional<OpenRing>, kMaxRings> rings_{};
    std::array<std::size_t, kMaxRings> ring_columns_{};
};

int bond_order_sum(const Graph& graph, std::size_t atom) noexcept {
    int sum = 0;
    for (const auto& bond : graph.bonds) {
        if (bond.first == atom || bond.second == atom) {
            sum += static_cast<int>(bond.order);
        }
    }
    return sum;
}

// Unbracketed atoms carry implicit hydrogens up to the lowest normal valence that fits.
void add_implicit_hydrogens(Graph& graph) {
    for (std::size_t i = 0; i < graph.atoms.size(); ++i) {
        auto& atom = graph.atoms[i];
        if (atom.bracket) {
            continue;
        }
        const int used = bond_order_sum(graph, i);
        for (const int valence : info(atom.element).valences) {
            if (valence != 0 && valence >= used) {
                atom.hydrogens = static_cast<std::uint8_t>(valence - used);
                break;
            }
        }
    }
}

}  // namespace

std::string_view symbol(Element element) noexcept {
    return info(element).symbol;
}

std::span<const std::size_t> Molecule::bonds_of(std::size_t atom) const noexcept {
    return std::span<const std::size_t>{adjacency_}.subspan(offsets_[atom],
                                                            offsets_[atom + 1] - offsets_[atom]);
}

std::size_t Molecule::neighbor(std::size_t bond, std::size_t atom) const noexcept {
    const auto& b = bonds_[bond];
    return b.first == atom ? b.second : b.first;
}

int Molecule::hydrogen_count(std::size_t atom) const noexcept {
    int count = atoms_[atom].hydrogens;
    for (const auto bond : bonds_of(atom)) {
        count += atoms_[neighbor(bond, atom)].element == Element::H ? 1 : 0;
    }
    return count;
}

int Molecule::valence(std::size_t atom) const noexcept {
    int total = atoms_[atom].hydrogens;
    for (const auto bond : bonds_of(atom)) {
        total += static_cast<int>(bonds_[bond].order);
    }
    return total;
}

core::Result<Molecule> parse_smiles(std::string_view smiles) {
    auto graph = Parser{smiles}.run();
    if (!graph.ok()) {
        return std::move(graph).error();
    }
    add_implicit_hydrogens(graph.value());

    Molecule molecule;
    molecule.atoms_ = std::move(graph.value().atoms);
    molecule.bonds_ = std::move(graph.value().bonds);
    const std::size_t n = molecule.atoms_.size();
    molecule.offsets_.assign(n + 1, 0);
    for (const auto& bond : molecule.bonds_) {
        ++molecule.offsets_[bond.first + 1];
        ++molecule.offsets_[bond.second + 1];
    }
    for (std::size_t i = 0; i < n; ++i) {
        molecule.offsets_[i + 1] += molecule.offsets_[i];
    }
    molecule.adjacency_.resize(molecule.offsets_[n]);
    std::vector<std::size_t> fill(molecule.offsets_.begin(), molecule.offsets_.end() - 1);
    for (std::size_t b = 0; b < molecule.bonds_.size(); ++b) {
        molecule.adjacency_[fill[molecule.bonds_[b].first]++] = b;
        molecule.adjacency_[fill[molecule.bonds_[b].second]++] = b;
    }
    return molecule;
}

std::string_view to_string(StructureIssue issue) noexcept {
    switch (issue) {
        case StructureIssue::Charged:
            return "charged";
        case StructureIssue::Radical:
            return "radical";
        case StructureIssue::UnusualValence:
            return "unusual valence";
    }
    return "?";
}

namespace {

// Neutral N at valence 5 with only single bonds: valence 5 needs a double bond (N(=O)=O).
bool is_five_sigma_nitrogen(const Molecule& molecule, std::size_t atom) {
    const auto& a = molecule.atoms()[atom];
    if (a.element != Element::N || a.charge != 0 || molecule.valence(atom) != 5) {
        return false;
    }
    const auto bonds = molecule.bonds_of(atom);
    return std::all_of(bonds.begin(), bonds.end(), [&](std::size_t bond) {
        return molecule.bonds()[bond].order == BondOrder::Single;
    });
}

}  // namespace

std::optional<StructureIssue> structure_issue(const Molecule& molecule) {
    const auto atoms = molecule.atoms();
    int net = 0;
    for (const auto& atom : atoms) {
        net += atom.charge;
    }
    if (net != 0) {
        return StructureIssue::Charged;
    }
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const int valence = molecule.valence(i);
        const auto normal = normal_valences(atoms[i]);
        if (!is_normal_valence(normal, valence) && valence < max_valence(normal)) {
            return StructureIssue::Radical;
        }
    }
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (molecule.valence(i) > max_valence(normal_valences(atoms[i])) ||
            is_five_sigma_nitrogen(molecule, i)) {
            return StructureIssue::UnusualValence;
        }
    }
    return std::nullopt;
}

}  // namespace arise::chem
