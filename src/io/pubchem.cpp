#include "arise/io/pubchem.hpp"

#include "arise/core/error.hpp"
#include "csv_table.hpp"
#include "smiles_text.hpp"
#include "text_lines.hpp"

#include <algorithm>
#include <optional>
#include <set>
#include <span>
#include <utility>

namespace arise::io {
namespace {

constexpr std::string_view kSearchBase =
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/";
constexpr std::string_view kSearchProperty = "/property/ConnectivitySMILES/CSV";
constexpr std::string_view kHeader = R"("CID","ConnectivitySMILES")";

struct Record {
    std::uint64_t cid;
    std::string smiles;
    std::size_t line;
};

core::Error parse_error(std::size_t line, std::string_view what) {
    return detail::error_at(core::ErrorCode::ParseError, line, what);
}

// A positive CID; up to 18 digits, so the value cannot overflow 64 bits.
std::optional<std::uint64_t> parse_cid(std::string_view text) noexcept {
    if (text.empty() || text.size() > 18) {
        return std::nullopt;
    }
    std::uint64_t value = 0;
    for (const char c : text) {
        if (c < '0' || c > '9') {
            return std::nullopt;
        }
        value = (value * 10U) + static_cast<std::uint64_t>(c - '0');
    }
    if (value == 0) {
        return std::nullopt;
    }
    return value;
}

core::Result<Record> parse_record(const detail::Line& line) {
    const auto text = detail::trim(line.text);
    const auto comma = text.find(',');
    if (comma == std::string_view::npos) {
        return parse_error(line.number, "expected CID,\"SMILES\"");
    }
    const auto cid_text = text.substr(0, comma);
    const auto cid = parse_cid(cid_text);
    if (!cid.has_value()) {
        return parse_error(line.number, "CID must be a positive whole number, got '" +
                                            std::string{cid_text} + "'");
    }
    auto field = text.substr(comma + 1);
    if (field.size() < 2 || field.front() != '"' || field.back() != '"') {
        return parse_error(line.number, "SMILES must be a quoted field");
    }
    field.remove_prefix(1);
    field.remove_suffix(1);
    if (field.empty()) {
        return parse_error(line.number, "empty SMILES for CID " + std::to_string(*cid));
    }
    if (const auto bad = detail::first_non_smiles_char(field)) {
        return parse_error(line.number, detail::describe_char(*bad) +
                                            " is not used in SMILES (CID " + std::to_string(*cid) +
                                            ")");
    }
    return Record{.cid = *cid, .smiles = std::string{field}, .line = line.number};
}

// '+' or '-' inside a bracket atom; outside brackets '-' is a single bond.
bool has_charged_atom(std::string_view smiles) noexcept {
    bool in_bracket = false;
    for (const char c : smiles) {
        if (c == '[') {
            in_bracket = true;
        } else if (c == ']') {
            in_bracket = false;
        } else if (in_bracket && (c == '+' || c == '-')) {
            return true;
        }
    }
    return false;
}

}  // namespace

std::string pubchem_formula_search_url(const core::Formula& formula) {
    return std::string{kSearchBase} + formula.hill() + std::string{kSearchProperty};
}

std::string_view to_string(RemovalReason reason) noexcept {
    switch (reason) {
        case RemovalReason::Mixture:
            return "mixture";
        case RemovalReason::Charged:
            return "charged";
        case RemovalReason::DuplicateConnectivity:
            return "duplicate structure";
    }
    return "?";
}

std::size_t PubChemImport::removed_count(RemovalReason reason) const noexcept {
    return static_cast<std::size_t>(
        std::count_if(removed.begin(), removed.end(),
                      [reason](const Removal& removal) { return removal.reason == reason; }));
}

core::Result<PubChemImport> import_pubchem_csv(std::string_view text) {
    const auto lines = detail::content_lines(text);
    if (lines.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "empty response: expected the header " + std::string{kHeader},
        };
    }
    if (detail::trim(lines.front().text) != kHeader) {
        return parse_error(lines.front().number,
                           "not a PubChem CID/ConnectivitySMILES CSV: expected the header " +
                               std::string{kHeader});
    }

    std::vector<Record> records;
    records.reserve(lines.size() - 1);
    for (const auto& line : std::span{lines}.subspan(1)) {
        auto record = parse_record(line);
        if (!record.ok()) {
            return std::move(record).error();
        }
        records.push_back(std::move(record).value());
    }
    // By CID, so the kept duplicate and every order below are independent of row order.
    std::stable_sort(records.begin(), records.end(),
                     [](const Record& a, const Record& b) { return a.cid < b.cid; });

    PubChemImport result;
    result.rows = records.size();
    std::set<std::string, std::less<>> seen;
    for (auto& record : records) {
        std::optional<RemovalReason> reason;
        if (record.smiles.find('.') != std::string::npos) {
            reason = RemovalReason::Mixture;
        } else if (has_charged_atom(record.smiles)) {
            reason = RemovalReason::Charged;
        } else if (!seen.insert(record.smiles).second) {
            reason = RemovalReason::DuplicateConnectivity;
        }
        if (reason.has_value()) {
            result.removed.push_back(
                Removal{.cid = record.cid, .smiles = std::move(record.smiles), .reason = *reason});
        } else {
            result.candidates.push_back(Candidate{
                .smiles = std::move(record.smiles),
                .name = "CID " + std::to_string(record.cid),
                .line = record.line,
            });
        }
    }
    return result;
}

}  // namespace arise::io
