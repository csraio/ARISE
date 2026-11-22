#pragma once

#include "prune_pipeline.hpp"

#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Development-only: the results the speed work of plan 3d must not change. For every
// committed PubChem search, each cleaned candidate's chem::nmr_equivalence (fewest and
// most); for the benchmark case's search, also its carbon fit in every scenario. The
// committed snapshot keeps one fingerprint per candidate, to stay small.
namespace arise::bench {

/// One recorded PubChem search: its formula and CSV text.
using Search = std::pair<std::string, std::string_view>;

/// Splits text of "## <formula>" lines, each followed by that search's CSV (the test
/// fixtures' kPubchemSearchesChunks), into searches in the order given; the views point into
/// `all`. Text before the first "## " line is ignored.
[[nodiscard]] std::vector<Search> split_searches(std::string_view all);

/// The provenance header of tests/golden/equivalence_snapshot.tsv (comment lines).
[[nodiscard]] std::string_view snapshot_header();

/// Every candidate's full result line, per search (in the order given): per cleaned
/// candidate (import order) its SMILES, fewest and most classes (space-separated); for the
/// search whose formula is the case's, one more field per scenario: tier, assumed overlaps
/// and reasons joined by " | ". ParseError if a search does not import, or as evidence_for.
struct SearchResults {
    std::string formula;
    std::vector<std::string> lines;
};
[[nodiscard]] core::Result<std::vector<SearchResults>> results(const std::vector<Search>& searches,
                                                               const BenchmarkCase& benchmark_case);

/// 64-bit FNV-1a of `text`, as 16 lowercase hex digits: a fingerprint of one result line.
[[nodiscard]] std::string fingerprint(std::string_view text);

/// The snapshot body: per search a "## <formula> <candidates>" line, then one fingerprint
/// per candidate's result line.
[[nodiscard]] std::string snapshot(const std::vector<SearchResults>& all);

}  // namespace arise::bench
