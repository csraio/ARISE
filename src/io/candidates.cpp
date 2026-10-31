#include "arise/io/candidates.hpp"

#include "arise/core/error.hpp"
#include "csv_table.hpp"
#include "smiles_text.hpp"
#include "text_lines.hpp"

namespace arise::io {
namespace {

constexpr std::string_view kBlank = " \t";

bool is_header_word(std::string_view token) noexcept {
    constexpr std::string_view kHeader = "smiles";
    if (token.size() != kHeader.size()) {
        return false;
    }
    for (std::size_t i = 0; i < token.size(); ++i) {
        const char lower =
            (token[i] >= 'A' && token[i] <= 'Z') ? static_cast<char>(token[i] + 32) : token[i];
        if (lower != kHeader[i]) {
            return false;
        }
    }
    return true;
}

core::Result<Candidate> parse_line(const detail::Line& line, bool first) {
    const auto content = detail::trim(line.text);
    const auto split = content.find_first_of(kBlank);
    const auto smiles = content.substr(0, split);
    const auto name =
        split == std::string_view::npos ? std::string_view{} : detail::trim(content.substr(split));

    if (const auto bad = detail::first_non_smiles_char(smiles)) {
        std::string what = detail::describe_char(*bad) + " is not used in SMILES";
        if (*bad == ',') {
            what += "; is this a CSV file? Separate SMILES and name with a space";
        }
        return detail::error_at(core::ErrorCode::ParseError, line.number, what);
    }
    if (first && is_header_word(smiles)) {
        return detail::error_at(
            core::ErrorCode::ParseError, line.number,
            "'" + std::string{smiles} + "' looks like a column header; SMILES files have none");
    }
    return Candidate{.smiles = std::string{smiles}, .name = std::string{name}, .line = line.number};
}

}  // namespace

core::Result<std::vector<Candidate>> parse_candidate_list(std::string_view text) {
    const auto lines = detail::content_lines(text);
    if (lines.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "no candidates: expected one SMILES per line",
        };
    }
    std::vector<Candidate> candidates;
    candidates.reserve(lines.size());
    for (const auto& line : lines) {
        auto candidate = parse_line(line, candidates.empty());
        if (!candidate.ok()) {
            return std::move(candidate).error();
        }
        candidates.push_back(std::move(candidate).value());
    }
    return candidates;
}

}  // namespace arise::io
