#include "text_lines.hpp"

namespace arise::io::detail {

std::string_view trim(std::string_view text) noexcept {
    constexpr std::string_view kBlank = " \t";
    const auto first = text.find_first_not_of(kBlank);
    if (first == std::string_view::npos) {
        return {};
    }
    text.remove_suffix(text.size() - 1 - text.find_last_not_of(kBlank));
    text.remove_prefix(first);
    return text;
}

std::vector<Line> content_lines(std::string_view text) {
    constexpr std::string_view kByteOrderMark = "\xEF\xBB\xBF";
    if (text.starts_with(kByteOrderMark)) {
        text.remove_prefix(kByteOrderMark.size());
    }
    std::vector<Line> lines;
    std::size_t number = 0;
    while (!text.empty()) {
        ++number;
        const auto newline = text.find('\n');
        std::string_view raw = text.substr(0, newline);
        text = newline == std::string_view::npos ? std::string_view{} : text.substr(newline + 1);
        if (!raw.empty() && raw.back() == '\r') {
            raw.remove_suffix(1);
        }
        const auto content = trim(raw);
        if (content.empty() || content.front() == '#') {
            continue;
        }
        lines.push_back(Line{.number = number, .text = raw});
    }
    return lines;
}

}  // namespace arise::io::detail
