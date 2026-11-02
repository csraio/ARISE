#pragma once

#include "arise/core/formula.hpp"
#include "arise/core/result.hpp"
#include "arise/io/candidates.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace arise::io {

/// The PubChem PUG REST address that lists every compound with exactly `formula`, as a
/// CSV of CID and ConnectivitySMILES. The engine only builds the address; the client
/// does the request.
[[nodiscard]] std::string pubchem_formula_search_url(const core::Formula& formula);

/// Why a PubChem record was not kept as a candidate.
enum class RemovalReason : std::uint8_t {
    Mixture,                // more than one component ('.')
    Charged,                // an ion: formal charges that do not sum to zero
    DuplicateConnectivity,  // same structure as a record with a lower CID
};

/// "mixture", "charged" or "duplicate structure".
[[nodiscard]] std::string_view to_string(RemovalReason reason) noexcept;

struct Removal {
    std::uint64_t cid = 0;
    std::string smiles;
    RemovalReason reason = RemovalReason::Mixture;

    friend bool operator==(const Removal&, const Removal&) = default;
};

/// The cleaned result of a PubChem formula search.
struct PubChemImport {
    /// One per distinct structure, by ascending CID. `name` is "CID <n>"; `line` is
    /// the CSV line of that record.
    std::vector<Candidate> candidates;
    /// Every record not kept, by ascending CID.
    std::vector<Removal> removed;
    /// Data rows read.
    std::size_t rows = 0;

    /// How many records were removed for `reason`.
    [[nodiscard]] std::size_t removed_count(RemovalReason reason) const noexcept;
};

/// Reads PubChem's CID,ConnectivitySMILES CSV and cleans it, removing, in this order:
/// mixtures, ions, then duplicate structures (keeping the lowest CID). Charge-separated
/// neutral molecules, whose formal charges sum to zero, are kept: PubChem writes nitro
/// groups, N-oxides, zwitterions and isocyanides that way ([N+](=O)[O-]).
/// Records are processed by ascending CID, so the result does not depend on PubChem's
/// row order. A header with no rows is a valid, empty result (no compound has that
/// formula). ParseError (naming the line) for anything that is not that CSV: a wrong
/// header, a bad CID, an unquoted or non-SMILES field.
[[nodiscard]] core::Result<PubChemImport> import_pubchem_csv(std::string_view text);

}  // namespace arise::io
