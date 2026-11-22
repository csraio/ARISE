#include "prune_pipeline.hpp"

#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"
#include "arise/io/number.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace arise::bench {

namespace {

std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    while (true) {
        const auto end = text.find(separator);
        parts.push_back(text.substr(0, end));
        if (end == std::string_view::npos) {
            return parts;
        }
        text = text.substr(end + 1);
    }
}

core::Error case_error(std::size_t line, const std::string& message) {
    return {
        .code = core::ErrorCode::ParseError,
        .message = "benchmark case, line " + std::to_string(line) + ": " + message,
    };
}

core::Error invalid(const std::string& message) {
    return {.code = core::ErrorCode::InvalidArgument, .message = "arise_bench: " + message};
}

std::optional<nmr1d::CarbonTypeSet> parse_types(std::string_view text) {
    constexpr std::array kTypes{
        nmr1d::CarbonType::C,
        nmr1d::CarbonType::CH,
        nmr1d::CarbonType::CH2,
        nmr1d::CarbonType::CH3,
    };
    nmr1d::CarbonTypeSet types;
    for (const auto name : split(text, ' ')) {
        bool known = false;
        for (const auto type : kTypes) {
            if (nmr1d::to_string(type) == name) {
                types = types | nmr1d::CarbonTypeSet{type};
                known = true;
            }
        }
        if (!known) {
            return std::nullopt;
        }
    }
    return types;
}

std::optional<long long> parse_whole(std::string_view text) {
    if (text.empty() || text.size() > 18 ||
        !std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; })) {
        return std::nullopt;
    }
    long long value = 0;
    for (const char c : text) {
        value = (value * 10) + (c - '0');
    }
    return value;
}

// The fields a case must have once each.
struct CaseFields {
    std::optional<std::string> formula;
    std::optional<long long> record;
    std::optional<std::string> structure;
};

// Reads one data line into `fields` or `c.peaks`.
std::optional<core::Error> read_case_line(std::size_t number, std::string_view line,
                                          CaseFields& fields, BenchmarkCase& c) {
    const auto parts = split(line, '\t');
    const auto kind = parts.front();
    const auto twice = [&] { return case_error(number, std::string{kind} + " listed twice"); };
    if (kind == "peak") {
        if (parts.size() != 3) {
            return case_error(number, "a peak needs a shift and its DEPT types");
        }
        const auto shift = io::parse_decimal(parts[1]);
        const auto types = parse_types(parts[2]);
        if (!shift || !std::isfinite(shift.value()) || !types) {
            return case_error(number, "bad peak '" + std::string{line} + "'");
        }
        c.peaks.push_back({.shift = core::Ppm{shift.value()}, .possible = *types});
        return std::nullopt;
    }
    if (parts.size() != 2 || parts[1].empty()) {
        return case_error(number, "expected '" + std::string{kind} + "' and one value");
    }
    if (kind == "formula") {
        if (fields.formula) {
            return twice();
        }
        fields.formula = std::string{parts[1]};
    } else if (kind == "record") {
        const auto record = parse_whole(parts[1]);
        if (!record) {
            return case_error(number, "bad record id '" + std::string{parts[1]} + "'");
        }
        if (fields.record) {
            return twice();
        }
        fields.record = record;
    } else if (kind == "structure") {
        if (fields.structure) {
            return twice();
        }
        fields.structure = std::string{parts[1]};
    } else {
        return case_error(number, "unknown line kind '" + std::string{kind} + "'");
    }
    return std::nullopt;
}

// The number of `most` carbon environments of `molecule` whose carbons are all sp.
int sp_environments(const chem::Molecule& molecule) {
    const auto classes = chem::nmr_equivalence(molecule).most;
    std::map<std::size_t, bool> all_sp;
    for (std::size_t atom = 0; atom < classes.size(); ++atom) {
        if (molecule.atoms()[atom].element != chem::Element::C) {
            continue;
        }
        const bool sp = chem::carbon_hybridization(molecule, atom) == chem::Hybridization::Sp;
        const auto [it, added] = all_sp.emplace(classes[atom], sp);
        if (!added) {
            it->second = it->second && sp;
        }
    }
    return static_cast<int>(std::count_if(all_sp.begin(), all_sp.end(),
                                          [](const auto& entry) { return entry.second; }));
}

}  // namespace

core::Result<BenchmarkCase> parse_benchmark_case(std::string_view text) {
    BenchmarkCase c;
    CaseFields fields;
    std::size_t number = 0;
    for (const auto line : split(text, '\n')) {
        ++number;
        if (line.empty() || line.front() == '#') {
            continue;
        }
        if (const auto error = read_case_line(number, line, fields, c)) {
            return *error;
        }
    }
    if (!fields.formula || !fields.record || !fields.structure || c.peaks.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "benchmark case: needs a formula, a record, a structure and peaks",
        };
    }
    c.formula = *fields.formula;
    c.record = *fields.record;
    c.structure = *fields.structure;
    return c;
}

std::string_view to_string(Scenario scenario) noexcept {
    switch (scenario) {
        case Scenario::Worst:
            return "worst";
        case Scenario::Typical:
            return "typical";
        case Scenario::NoDept:
            return "no-dept";
    }
    return "";
}

core::Result<Evidence> evidence_for(Scenario scenario, const BenchmarkCase& c) {
    Evidence evidence;
    evidence.peaks = c.peaks;
    if (scenario == Scenario::Worst) {
        return evidence;
    }
    const auto truth = chem::parse_smiles(c.structure);
    if (!truth) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "benchmark case: the structure does not read: " + truth.error().message,
        };
    }
    if (scenario == Scenario::Typical) {
        for (const auto group : chem::kAllGroups) {
            if (chem::count_group(truth.value(), group) > 0) {
                evidence.groups.present.push_back(group);
            }
        }
        const int exchangeable = chem::exchangeable_hydrogens(truth.value());
        evidence.exchangeable =
            filter::ConfirmedExchangeable{.min = exchangeable, .max = exchangeable};
        return evidence;
    }
    for (auto& peak : evidence.peaks) {
        peak.possible = nmr1d::CarbonTypeSet::all();
    }
    const auto windows = filter::measured_hybridization_windows();
    if (!windows) {
        return windows.error();
    }
    const int sp = sp_environments(truth.value());
    const auto suggestion = filter::suggest_hybridization(evidence.peaks, windows.value(), sp);
    if (!suggestion) {
        return suggestion.error();
    }
    filter::ConfirmedCarbonCounts counts;
    if (sp > 0) {
        counts.sp = sp;
    }
    counts.sp2 = suggestion.value().sp2.min;
    counts.sp3 = suggestion.value().sp3.min;
    evidence.counts = counts;
    return evidence;
}

core::Result<io::PubChemImport> clean(std::string_view pubchem_csv) {
    return io::import_pubchem_csv(pubchem_csv);
}

Parsed parse(const io::PubChemImport& imported) {
    Parsed parsed;
    parsed.molecules.reserve(imported.candidates.size());
    for (const auto& candidate : imported.candidates) {
        auto molecule = chem::parse_smiles(candidate.smiles);
        if (molecule) {
            parsed.molecules.push_back(std::move(molecule).value());
            parsed.smiles.push_back(candidate.smiles);
        } else {
            ++parsed.unreadable;
        }
    }
    return parsed;
}

core::Result<std::vector<filter::CarbonFit>> prune(std::span<const chem::Molecule> molecules,
                                                   const Evidence& evidence) {
    std::vector<filter::CarbonFit> fits;
    fits.reserve(molecules.size());
    for (const auto& molecule : molecules) {
        auto groups = filter::check_confirmed_groups(molecule, evidence.groups);
        if (!groups) {
            return groups.error();
        }
        if (groups.value().tier == filter::Tier::Active && evidence.exchangeable) {
            groups = filter::check_confirmed_exchangeable(molecule, *evidence.exchangeable);
            if (!groups) {
                return groups.error();
            }
        }
        if (groups.value().tier != filter::Tier::Active) {
            filter::CarbonFit eliminated;
            eliminated.verdict = std::move(groups).value();
            fits.push_back(std::move(eliminated));
            continue;
        }
        auto fit = evidence.counts
                       ? filter::check_carbon_fit(molecule, evidence.peaks, *evidence.counts)
                       : filter::check_carbon_fit(molecule, evidence.peaks);
        if (!fit) {
            return fit.error();
        }
        fits.push_back(std::move(fit).value());
    }
    return fits;
}

std::vector<std::size_t> rank(std::span<const filter::CarbonFit> fits) {
    return filter::order_by_fit(fits);
}

core::Result<Summary> summarize(std::vector<double> samples_ms) {
    if (samples_ms.empty()) {
        return invalid("no timing samples");
    }
    if (std::any_of(samples_ms.begin(), samples_ms.end(),
                    [](double s) { return !std::isfinite(s) || s < 0; })) {
        return invalid("a timing sample is negative or not finite");
    }
    std::sort(samples_ms.begin(), samples_ms.end());
    const std::size_t n = samples_ms.size();
    const double median =
        n % 2 == 1 ? samples_ms[n / 2] : (samples_ms[(n / 2) - 1] + samples_ms[n / 2]) / 2.0;
    return Summary{.min = samples_ms.front(), .median = median, .max = samples_ms.back()};
}

core::Result<Options> parse_options(std::span<const std::string_view> args) {
    Options options;
    for (std::size_t i = 0; i < args.size(); i += 2) {
        if (i + 1 >= args.size()) {
            return invalid("'" + std::string{args[i]} + "' needs a value");
        }
        const auto value = args[i + 1];
        if (args[i] == "--runs") {
            const auto runs = parse_whole(value);
            if (!runs || *runs < 1) {
                return invalid("--runs needs a whole number of at least 1, got '" +
                               std::string{value} + "'");
            }
            options.runs = static_cast<std::size_t>(*runs);
        } else if (args[i] == "--budget-ms") {
            const auto budget = io::parse_decimal(value);
            if (!budget || !std::isfinite(budget.value()) || budget.value() <= 0) {
                return invalid("--budget-ms needs a positive number, got '" + std::string{value} +
                               "'");
            }
            options.budget_ms = budget.value();
        } else {
            return invalid("unknown option '" + std::string{args[i]} + "'");
        }
    }
    return options;
}

int exit_code(double median_ms, std::optional<double> budget_ms) noexcept {
    return budget_ms && median_ms > *budget_ms ? 1 : 0;
}

}  // namespace arise::bench
