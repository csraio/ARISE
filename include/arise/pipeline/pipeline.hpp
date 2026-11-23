#pragma once
// The core loop of the web client (plan "web client", 2026-11-22), as stateless calls:
// every call takes all its input as text or plain values and returns a value, so the same
// input gives the same output, natively and in WebAssembly. The page keeps the state.
//
//   formula -> search_url -> (the page fetches PubChem) -> rank
//   spectra -> read_spectra / suggest -> (the chemist confirms groups) -> rank

#include "arise/chem/perception.hpp"
#include "arise/core/result.hpp"
#include "arise/filter/groups.hpp"
#include "arise/io/pubchem.hpp"
#include "arise/nmr1d/dept.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace arise::pipeline {

/// How a spectrum's text is laid out: ARISE's peak-list CSV (an uploaded file, or the
/// page's table written out as CSV; io::parse_carbon_csv, io::parse_dept_csv) or a list
/// pasted by the chemist (io::parse_pasted_carbon, io::parse_pasted_dept).
enum class TextFormat : std::uint8_t { Csv, Pasted };

/// "csv" or "pasted".
[[nodiscard]] std::string_view to_string(TextFormat format) noexcept;

/// The format whose to_string is `name`, if any.
[[nodiscard]] std::optional<TextFormat> text_format_named(std::string_view name);

struct SpectrumText {
    TextFormat format = TextFormat::Pasted;
    std::string text;
};

/// The spectra as the chemist entered them: a 13C list, and DEPT-135 and DEPT-90 when
/// recorded. An absent DEPT spectrum differs from one without peaks.
struct SpectraText {
    SpectrumText carbon;
    std::optional<SpectrumText> dept135;
    std::optional<SpectrumText> dept90;
};

/// The PubChem address listing every compound with `formula` (io::pubchem_formula_search_url).
/// ParseError when the formula does not read.
[[nodiscard]] core::Result<std::string> search_url(std::string_view formula);

/// Reads the spectra and classifies each 13C peak by its DEPT phases
/// (nmr1d::classify_carbons; every type is possible without DEPT). Errors keep their code
/// and name the spectrum: "13C spectrum: line 3: ...", "DEPT-135 spectrum: ...",
/// "DEPT-90 spectrum: ...".
[[nodiscard]] core::Result<std::vector<nmr1d::CarbonAssignment>> read_spectra(
    const SpectraText& spectra);

/// The groups the spectra suggest from ARISE's built-in measured windows
/// (filter::suggest_groups with filter::measured_hallmark_windows). Suggestions are never
/// eliminations; the chemist confirms.
[[nodiscard]] core::Result<std::vector<filter::GroupSuggestion>> suggest(
    const SpectraText& spectra, filter::SuggestOptions options);

/// The group whose chem::to_string is `name` ("ketone"), if any.
[[nodiscard]] std::optional<chem::Group> group_named(std::string_view name);

/// One candidate's place in the ranking.
struct RankedCandidate {
    std::uint64_t cid = 0;
    std::string smiles;
    /// Tier and reasons: Eliminated when the SMILES does not read or a confirmed group is
    /// violated (filter::check_confirmed_groups); otherwise the carbon fit's verdict
    /// (filter::check_carbon_fit), at least BackBurner when a group confirmed only as likely
    /// (from 1H) is missing, that group's reasons first.
    filter::Verdict verdict;
    /// From the carbon fit; zero for a candidate eliminated before it.
    int assumed_overlaps = 0;
    int environments_fewest = 0;
    int environments_most = 0;
    std::vector<std::string> notes;

    friend bool operator==(const RankedCandidate&, const RankedCandidate&) = default;
};

/// The cleaned PubChem list, ranked against the spectra and the confirmed groups.
struct Ranking {
    /// The formula searched, in Hill order.
    std::string formula;
    /// PubChem rows read, and the records cleaning removed (io::import_pubchem_csv).
    std::size_t rows = 0;
    std::vector<io::Removal> removed;
    /// Every cleaned candidate, best first (filter::order_by_fit: tier, then assumed
    /// overlaps, then number of reasons; ties keep ascending CID).
    std::vector<RankedCandidate> candidates;

    friend bool operator==(const Ranking&, const Ranking&) = default;
};

/// Cleans `pubchem_csv` (PubChem's CID,ConnectivitySMILES answer for `formula`) and ranks
/// every candidate: confirmed groups first, then the carbon fit against the 13C peaks and
/// their DEPT types. Nothing is deleted; every candidate keeps its reasons.
/// Errors: ParseError for a formula, CSV or spectrum that does not read; InvalidArgument
/// for a group confirmed both present and absent.
[[nodiscard]] core::Result<Ranking> rank(std::string_view formula, std::string_view pubchem_csv,
                                         const SpectraText& spectra,
                                         const filter::ConfirmedGroups& confirmed,
                                         filter::SuggestOptions options);

}  // namespace arise::pipeline
