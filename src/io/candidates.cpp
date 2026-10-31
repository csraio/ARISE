#include "arise/io/candidates.hpp"

#include "arise/core/error.hpp"
#include "csv_table.hpp"
#include "text_lines.hpp"

namespace arise::io {
namespace {

constexpr std::string_view kBlank = " \t";

bool is_smiles_char(char c) noexcept {
    if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')) {
        return true;
    }
    constexpr std::string_view kSymbols = "()[]=#$:/\\.+-@%*";
    return kSymbols.find(c) != std::string_view::npos;
}

// How to name a rejected character in a message.
std::string describe(char c) {
    const auto byte = static_cast<unsigned char>(c);
    if (byte >= 0x80) {
        return "a non-ASCII byte";
    }
    if (byte < 0x20 || byte == 0x7F) {
        constexpr std::string_view kHexDigits = "0123456789ABCDEF";
        return std::string{"control character 0x"} + kHexDigits[byte >> 4U] +
               kHexDigits[byte & 0xFU];
    }
    return "'" + std::string(1, c) + "'";
}

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

    for (const char c : smiles) {
        if (!is_smiles_char(c)) {
            std::string what = describe(c) + " is not used in SMILES";
            if (c == ',') {
                what += "; is this a CSV file? Separate SMILES and name with a space";
            }
            return detail::error_at(core::ErrorCode::ParseError, line.number, what);
        }
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
