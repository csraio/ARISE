#include "arise/nmr1d/solvent.hpp"

#include "arise/core/error.hpp"
#include "arise/io/number.hpp"
#include "arise/io/peak_csv.hpp"
#include "solvent_peaks_data.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>
#include <string>
#include <utility>

namespace arise::nmr1d {

namespace {

// Columns of data/solvent_peaks.tsv.
constexpr std::size_t kColumns = 5;
constexpr std::size_t kCidColumn = 1;
constexpr std::size_t kResidualColumn = 2;
constexpr std::size_t kWaterColumn = 3;
constexpr std::size_t kExchangeColumn = 4;

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

core::Error parse_error(std::size_t line, const std::string& message) {
    return {
        .code = core::ErrorCode::ParseError,
        .message = "solvent peaks, line " + std::to_string(line) + ": " + message,
    };
}

std::optional<Solvent> solvent_named(std::string_view name) {
    for (const auto solvent : kAllSolvents) {
        if (to_string(solvent) == name) {
            return solvent;
        }
    }
    return std::nullopt;
}

std::string solvent_names() {
    std::string names;
    for (const auto solvent : kAllSolvents) {
        names += names.empty() ? "" : ", ";
        names += to_string(solvent);
    }
    return names;
}

std::optional<std::size_t> parse_positive_count(std::string_view text) {
    if (text.empty() ||
        !std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; })) {
        return std::nullopt;
    }
    std::size_t value = 0;
    for (const char c : text) {
        value = (value * 10) + static_cast<std::size_t>(c - '0');
    }
    return value == 0 ? std::nullopt : std::optional{value};
}

// A comma-separated, ascending list of shifts in the 1H window.
core::Result<std::vector<core::Ppm>> parse_shifts(std::size_t line, std::string_view text) {
    std::vector<core::Ppm> shifts;
    for (const auto item : split(text, ',')) {
        const auto value = io::parse_decimal(item);
        if (!value || !std::isfinite(value.value())) {
            return parse_error(line, "not a number: '" + std::string{item} + "'");
        }
        const core::Ppm shift{value.value()};
        if (shift < io::kMinProtonShift || shift > io::kMaxProtonShift) {
            return parse_error(
                line, "shift " + std::string{item} + " is outside the 1H window [-10, 25] ppm");
        }
        if (!shifts.empty() && shift < shifts.back()) {
            return parse_error(line, "shift " + std::string{item} + " breaks ascending order");
        }
        shifts.push_back(shift);
    }
    return shifts;
}

core::Result<SolventPeaks> parse_row(std::size_t line, std::string_view content) {
    const auto fields = split(content, '\t');
    if (fields.size() != kColumns) {
        return parse_error(line, "expected 5 fields, got " + std::to_string(fields.size()));
    }
    const auto solvent = solvent_named(fields[0]);
    if (!solvent) {
        return parse_error(line, "unknown solvent '" + std::string{fields[0]} +
                                     "' (known: " + solvent_names() + ")");
    }
    const auto cid = parse_positive_count(fields[kCidColumn]);
    if (!cid) {
        return parse_error(line, "PubChem CID must be a positive whole number, got '" +
                                     std::string{fields[kCidColumn]} + "'");
    }
    if (fields[kResidualColumn] == "-") {
        return parse_error(line, "no residual shift");
    }
    auto residual = parse_shifts(line, fields[kResidualColumn]);
    if (!residual) {
        return std::move(residual).error();
    }
    std::vector<core::Ppm> water;
    if (fields[kWaterColumn] != "-") {
        auto parsed = parse_shifts(line, fields[kWaterColumn]);
        if (!parsed) {
            return std::move(parsed).error();
        }
        water = std::move(parsed).value();
    }
    const auto exchanges = fields[kExchangeColumn];
    if (exchanges != "yes" && exchanges != "no") {
        return parse_error(
            line, "exchanges labile H must be 'yes' or 'no', got '" + std::string{exchanges} + "'");
    }
    return SolventPeaks{
        .solvent = *solvent,
        .pubchem_cid = *cid,
        .residual = std::move(residual).value(),
        .water = std::move(water),
        .exchanges_labile_h = exchanges == "yes",
    };
}

}  // namespace

std::string_view to_string(Solvent solvent) noexcept {
    switch (solvent) {
        case Solvent::ThfD8:
            return "THF-d8";
        case Solvent::CD2Cl2:
            return "CD2Cl2";
        case Solvent::CDCl3:
            return "CDCl3";
        case Solvent::TolueneD8:
            return "toluene-d8";
        case Solvent::C6D6:
            return "C6D6";
        case Solvent::AcetoneD6:
            return "(CD3)2CO";
        case Solvent::DmsoD6:
            return "(CD3)2SO";
        case Solvent::CD3CN:
            return "CD3CN";
        case Solvent::TfeD3:
            return "TFE-d3";
        case Solvent::CD3OD:
            return "CD3OD";
        case Solvent::D2O:
            return "D2O";
    }
    return "unknown";
}

core::Result<std::vector<SolventPeaks>> parse_solvent_peaks(std::string_view text) {
    std::vector<SolventPeaks> table;
    std::set<Solvent> seen;
    std::size_t line = 0;
    for (auto content : split(text, '\n')) {
        ++line;
        if (content.ends_with('\r')) {
            content.remove_suffix(1);
        }
        if (content.empty() || content.front() == '#') {
            continue;
        }
        auto row = parse_row(line, content);
        if (!row) {
            return std::move(row).error();
        }
        if (!seen.insert(row.value().solvent).second) {
            return parse_error(line,
                               std::string{to_string(row.value().solvent)} + " is listed twice");
        }
        table.push_back(std::move(row).value());
    }
    if (table.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "solvent peaks: the table has no rows",
        };
    }
    return table;
}

core::Result<SolventPeaks> solvent_peaks(Solvent solvent) {
    auto table = parse_solvent_peaks(detail::kSolventPeaks);
    if (!table) {
        return std::move(table).error();
    }
    for (auto& row : table.value()) {
        if (row.solvent == solvent) {
            return std::move(row);
        }
    }
    return core::Error{
        .code = core::ErrorCode::NotFound,
        .message = "solvent peaks: no row for " + std::string{to_string(solvent)},
    };
}

}  // namespace arise::nmr1d
