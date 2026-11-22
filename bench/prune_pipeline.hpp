#pragma once

#include "arise/chem/molecule.hpp"
#include "arise/core/result.hpp"
#include "arise/filter/carbon_fit.hpp"
#include "arise/filter/exchangeable.hpp"
#include "arise/filter/groups.hpp"
#include "arise/filter/hybridization.hpp"
#include "arise/io/pubchem.hpp"
#include "arise/nmr1d/dept.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

// Development-only (never part of the engine): the cleaning and pruning chain a chemist's
// run takes (plan 3d), split into stages so bench/arise_bench can time each. It becomes the
// pipeline module later.
namespace arise::bench {

/// The benchmark case (tests/data/reference/benchmark_spectrum.tsv): a real nmrshiftdb2
/// record whose formula's PubChem search has 10,000+ rows.
struct BenchmarkCase {
    std::string formula;
    long long record = 0;
    std::string structure;                       // the true candidate's PubChem SMILES
    std::vector<nmr1d::CarbonAssignment> peaks;  // with the DEPT types of their carbons
};

/// Reads benchmark_spectrum.tsv: '#' comments, then tab-separated lines "formula F",
/// "record N", "structure SMILES" (once each) and "peak SHIFT TYPES" (TYPES: C, CH, CH2,
/// CH3 separated by spaces). ParseError, citing the line, for an unknown line kind, a
/// missing or repeated field, a bad number or type, or no peaks.
[[nodiscard]] core::Result<BenchmarkCase> parse_benchmark_case(std::string_view text);

/// What the chemist has confirmed, per benchmark scenario.
enum class Scenario : std::uint8_t {
    Worst,    // DEPT types, nothing else confirmed: every candidate reaches the carbon fit
              // (named for the evidence; the NoDept chain does more work per candidate)
    Typical,  // DEPT types, the true structure's groups confirmed present, exact OH/NH count
    NoDept,   // peaks without DEPT types, the accepted sp2/sp3 suggestions (typed sp count)
};

/// "worst", "typical", "no-dept".
[[nodiscard]] std::string_view to_string(Scenario scenario) noexcept;

/// The chemist's evidence for one scenario.
struct Evidence {
    std::vector<nmr1d::CarbonAssignment> peaks;
    filter::ConfirmedGroups groups;
    std::optional<filter::ConfirmedExchangeable> exchangeable;
    std::optional<filter::ConfirmedCarbonCounts> counts;
};

/// The evidence a chemist would give for `scenario`, derived from the case's true structure
/// with ARISE's own functions (inputs to time, not expected values). NoDept accepts
/// filter::suggest_hybridization's minimums with the measured windows, typing the true
/// structure's number of sp carbon environments. ParseError if the structure does not read.
[[nodiscard]] core::Result<Evidence> evidence_for(Scenario scenario, const BenchmarkCase& c);

/// Stage 1: io::import_pubchem_csv.
[[nodiscard]] core::Result<io::PubChemImport> clean(std::string_view pubchem_csv);

/// Stage 2: the candidates read as molecules, in import order; a candidate the reader
/// rejects is left out and counted in `unreadable`.
struct Parsed {
    std::vector<chem::Molecule> molecules;
    std::vector<std::string> smiles;  // of each molecule
    std::size_t unreadable = 0;
};
[[nodiscard]] Parsed parse(const io::PubChemImport& imported);

/// Stage 3: hard checks first (confirmed groups, then the confirmed exchangeable count);
/// a candidate they eliminate gets that verdict, and the carbon fit runs only on the rest
/// (check_carbon_fit with the confirmed counts when given). One fit per molecule.
[[nodiscard]] core::Result<std::vector<filter::CarbonFit>> prune(
    std::span<const chem::Molecule> molecules, const Evidence& evidence);

/// Stage 4: filter::order_by_fit.
[[nodiscard]] std::vector<std::size_t> rank(std::span<const filter::CarbonFit> fits);

/// Timing statistics of repeated runs, in milliseconds.
struct Summary {
    double min = 0;
    double median = 0;  // the mean of the middle two for an even count
    double max = 0;
};

/// InvalidArgument for no samples or a sample that is negative or not finite.
[[nodiscard]] core::Result<Summary> summarize(std::vector<double> samples_ms);

/// arise_bench's options: "--runs N" (N >= 1, default 15) and "--budget-ms X" (X > 0;
/// none by default). InvalidArgument for anything else.
struct Options {
    std::size_t runs = 15;
    std::optional<double> budget_ms;
};
[[nodiscard]] core::Result<Options> parse_options(std::span<const std::string_view> args);

/// 0 when there is no budget or `median_ms` is within it, else 1.
[[nodiscard]] int exit_code(double median_ms, std::optional<double> budget_ms) noexcept;

}  // namespace arise::bench
