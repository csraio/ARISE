#include "arise/pipeline/pipeline.hpp"

#include "arise/chem/molecule.hpp"
#include "arise/core/error.hpp"
#include "arise/core/formula.hpp"
#include "arise/filter/carbon_fit.hpp"
#include "arise/io/carbon_csv.hpp"
#include "arise/io/pasted_peaks.hpp"

#include <span>
#include <utility>

namespace arise::pipeline {
namespace {

// `error` with its message prefixed by "<what>: ".
core::Error named(std::string_view what, core::Error error) {
    error.message = std::string{what} + ": " + error.message;
    return error;
}

core::Result<std::vector<io::CarbonPeak>> read_carbon(const SpectrumText& spectrum) {
    auto peaks = spectrum.format == TextFormat::Csv ? io::parse_carbon_csv(spectrum.text)
                                                    : io::parse_pasted_carbon(spectrum.text);
    if (!peaks.ok()) {
        return named("13C spectrum", std::move(peaks).error());
    }
    return peaks;
}

core::Result<std::optional<std::vector<io::DeptPeak>>> read_dept(
    std::string_view name, const std::optional<SpectrumText>& spectrum) {
    if (!spectrum.has_value()) {
        return std::optional<std::vector<io::DeptPeak>>{};
    }
    auto peaks = spectrum->format == TextFormat::Csv ? io::parse_dept_csv(spectrum->text)
                                                     : io::parse_pasted_dept(spectrum->text);
    if (!peaks.ok()) {
        return named(name, std::move(peaks).error());
    }
    return std::optional{std::move(peaks).value()};
}

// The CID of a cleaned PubChem candidate, whose name is "CID <n>" (io::import_pubchem_csv).
core::Result<std::uint64_t> cid_of(const io::Candidate& candidate) {
    constexpr std::string_view kPrefix = "CID ";
    const std::string_view name = candidate.name;
    std::uint64_t cid = 0;
    if (!name.starts_with(kPrefix) || name.size() == kPrefix.size()) {
        return core::Error{
            .code = core::ErrorCode::Internal,
            .message = "PubChem candidate without a CID: '" + candidate.name + "'",
        };
    }
    for (const char c : name.substr(kPrefix.size())) {
        if (c < '0' || c > '9') {
            return core::Error{
                .code = core::ErrorCode::Internal,
                .message = "PubChem candidate without a CID: '" + candidate.name + "'",
            };
        }
        cid = (cid * 10) + static_cast<std::uint64_t>(c - '0');
    }
    return cid;
}

// The candidate's fit: Eliminated when its SMILES does not read or it violates a confirmed
// group; otherwise the carbon fit.
core::Result<filter::CarbonFit> fit_of(const io::Candidate& candidate,
                                       std::span<const nmr1d::CarbonAssignment> peaks,
                                       const filter::ConfirmedGroups& confirmed,
                                       filter::SuggestOptions options) {
    filter::CarbonFit fit;
    auto molecule = chem::parse_smiles(candidate.smiles);
    if (!molecule.ok()) {
        fit.verdict.tier = filter::Tier::Eliminated;
        fit.verdict.reasons.push_back("the SMILES does not read: " + molecule.error().message);
        return fit;
    }
    auto groups = filter::check_confirmed_groups(molecule.value(), confirmed);
    if (!groups.ok()) {
        return std::move(groups).error();
    }
    if (groups.value().tier == filter::Tier::Eliminated) {
        fit.verdict = std::move(groups).value();
        return fit;
    }
    return filter::check_carbon_fit(molecule.value(), peaks, options);
}

}  // namespace

core::Result<std::string> search_url(std::string_view formula) {
    const auto parsed = core::Formula::parse(formula);
    if (!parsed.ok()) {
        return parsed.error();
    }
    return io::pubchem_formula_search_url(parsed.value());
}

core::Result<std::vector<nmr1d::CarbonAssignment>> read_spectra(const SpectraText& spectra) {
    const auto carbon = read_carbon(spectra.carbon);
    if (!carbon.ok()) {
        return carbon.error();
    }
    const auto dept135 = read_dept("DEPT-135 spectrum", spectra.dept135);
    if (!dept135.ok()) {
        return dept135.error();
    }
    const auto dept90 = read_dept("DEPT-90 spectrum", spectra.dept90);
    if (!dept90.ok()) {
        return dept90.error();
    }
    nmr1d::DeptSpectra dept{.carbon = carbon.value()};
    if (dept135.value().has_value()) {
        dept.dept135 = std::span<const io::DeptPeak>{*dept135.value()};
    }
    if (dept90.value().has_value()) {
        dept.dept90 = std::span<const io::DeptPeak>{*dept90.value()};
    }
    return nmr1d::classify_carbons(dept);
}

core::Result<std::vector<filter::GroupSuggestion>> suggest(const SpectraText& spectra,
                                                           filter::SuggestOptions options) {
    const auto peaks = read_spectra(spectra);
    if (!peaks.ok()) {
        return peaks.error();
    }
    const auto windows = filter::measured_hallmark_windows();
    if (!windows.ok()) {
        return windows.error();
    }
    return filter::suggest_groups(peaks.value(), windows.value(), options);
}

std::optional<chem::Group> group_named(std::string_view name) {
    for (const auto group : chem::kAllGroups) {
        if (chem::to_string(group) == name) {
            return group;
        }
    }
    return std::nullopt;
}

core::Result<Ranking> rank(std::string_view formula, std::string_view pubchem_csv,
                           const SpectraText& spectra, const filter::ConfirmedGroups& confirmed,
                           filter::SuggestOptions options) {
    const auto parsed = core::Formula::parse(formula);
    if (!parsed.ok()) {
        return parsed.error();
    }
    const auto peaks = read_spectra(spectra);
    if (!peaks.ok()) {
        return peaks.error();
    }
    auto imported = io::import_pubchem_csv(pubchem_csv);
    if (!imported.ok()) {
        return std::move(imported).error();
    }
    const auto& candidates = imported.value().candidates;
    std::vector<filter::CarbonFit> fits;
    fits.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        auto fit = fit_of(candidate, peaks.value(), confirmed, options);
        if (!fit.ok()) {
            return std::move(fit).error();
        }
        fits.push_back(std::move(fit).value());
    }
    Ranking ranking{
        .formula = parsed.value().hill(),
        .rows = imported.value().rows,
        .removed = std::move(imported.value().removed),
        .candidates = {},
    };
    ranking.candidates.reserve(candidates.size());
    for (const auto index : filter::order_by_fit(fits)) {
        const auto cid = cid_of(candidates[index]);
        if (!cid.ok()) {
            return cid.error();
        }
        auto& fit = fits[index];
        ranking.candidates.push_back({
            .cid = cid.value(),
            .smiles = candidates[index].smiles,
            .verdict = std::move(fit.verdict),
            .assumed_overlaps = fit.assumed_overlaps,
            .environments_fewest = fit.environments_fewest,
            .environments_most = fit.environments_most,
            .notes = std::move(fit.notes),
        });
    }
    return ranking;
}

}  // namespace arise::pipeline
