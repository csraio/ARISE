#include "arise/core/formula.hpp"

#include "arise/core/error.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace arise::core {
namespace {

struct ElementInfo {
    std::string_view symbol;
    int valence;  // lowest common valence, for degrees of unsaturation
};

constexpr std::array kElements{
    ElementInfo{.symbol = "B", .valence = 3},  ElementInfo{.symbol = "Br", .valence = 1},
    ElementInfo{.symbol = "C", .valence = 4},  ElementInfo{.symbol = "Cl", .valence = 1},
    ElementInfo{.symbol = "F", .valence = 1},  ElementInfo{.symbol = "H", .valence = 1},
    ElementInfo{.symbol = "I", .valence = 1},  ElementInfo{.symbol = "N", .valence = 3},
    ElementInfo{.symbol = "O", .valence = 2},  ElementInfo{.symbol = "P", .valence = 3},
    ElementInfo{.symbol = "S", .valence = 2},  ElementInfo{.symbol = "Se", .valence = 2},
    ElementInfo{.symbol = "Si", .valence = 4},
};
constexpr int kMaxCount = 999;

const ElementInfo* find_element(std::string_view symbol) noexcept {
    for (const auto& element : kElements) {
        if (element.symbol == symbol) {
            return &element;
        }
    }
    return nullptr;
}

bool is_upper(char c) noexcept {
    return c >= 'A' && c <= 'Z';
}
bool is_lower(char c) noexcept {
    return c >= 'a' && c <= 'z';
}
bool is_digit(char c) noexcept {
    return c >= '0' && c <= '9';
}

Error parse_error(std::string message) {
    return Error{.code = ErrorCode::ParseError, .message = std::move(message)};
}

std::string_view trim(std::string_view text) noexcept {
    constexpr std::string_view kSpace = " \t\r\n";
    const auto first = text.find_first_not_of(kSpace);
    if (first == std::string_view::npos) {
        return {};
    }
    text.remove_suffix(text.size() - 1 - text.find_last_not_of(kSpace));
    text.remove_prefix(first);
    return text;
}

// Error for an unknown symbol, with a hint when it reads as two known one-letter
// elements whose second letter lost its capital ("Co" for "CO").
Error unknown_element(std::string_view symbol) {
    std::string message = "unsupported element '" + std::string{symbol} + "'";
    if (symbol.size() == 2) {
        const std::string first(1, symbol[0]);
        const std::string second(1, static_cast<char>(symbol[1] - 'a' + 'A'));
        if (find_element(first) != nullptr && find_element(second) != nullptr) {
            message +=
                "; did you mean '" + first + second + "'? Element symbols are case-sensitive";
        }
    }
    return parse_error(std::move(message));
}

void append_term(std::string& out, std::string_view symbol, int count) {
    out += symbol;
    if (count != 1) {
        out += std::to_string(count);
    }
}

// The element symbol starting at `pos` (one capital, optional lowercase letter).
Result<std::string_view> read_symbol(std::string_view text, std::size_t pos) {
    const char c = text[pos];
    if (!is_upper(c)) {
        if (is_lower(c) && pos == 0) {
            return parse_error(
                "formula must start with a capital letter (element symbols "
                "are case-sensitive), got '" +
                std::string{text} + "'");
        }
        if (is_digit(c)) {
            return parse_error("count '" + std::string(1, c) + "' must follow an element symbol");
        }
        return parse_error("unexpected character '" + std::string(1, c) + "' in '" +
                           std::string{text} + "'");
    }
    const std::size_t length = (pos + 1 < text.size() && is_lower(text[pos + 1])) ? 2 : 1;
    const auto symbol = text.substr(pos, length);
    if (find_element(symbol) == nullptr) {
        return unknown_element(symbol);
    }
    return symbol;
}

// The count written as `digits` after `symbol`; 1 when there are none.
Result<int> read_count(std::string_view digits, std::string_view symbol) {
    if (digits.empty()) {
        return 1;
    }
    if (digits.front() == '0') {
        return parse_error("count '" + std::string{digits} + "' after " + std::string{symbol} +
                           " must be a positive number");
    }
    if (digits.size() > 3) {
        return Error{
            .code = ErrorCode::OutOfRange,
            .message = "count " + std::string{digits} + " after " + std::string{symbol} +
                       " exceeds " + std::to_string(kMaxCount),
        };
    }
    int count = 0;
    for (const char digit : digits) {
        count = (count * 10) + (digit - '0');
    }
    return count;
}

}  // namespace

Result<Formula> Formula::parse(std::string_view text) {
    text = trim(text);
    if (text.empty()) {
        return parse_error("empty formula");
    }
    Formula formula;
    std::size_t pos = 0;
    while (pos < text.size()) {
        auto symbol = read_symbol(text, pos);
        if (!symbol.ok()) {
            return std::move(symbol).error();
        }
        pos += symbol.value().size();
        std::size_t end = pos;
        while (end < text.size() && is_digit(text[end])) {
            ++end;
        }
        auto count = read_count(text.substr(pos, end - pos), symbol.value());
        if (!count.ok()) {
            return std::move(count).error();
        }
        pos = end;
        auto& total = formula.counts_[std::string{symbol.value()}];
        total += count.value();
        if (total > kMaxCount) {
            return Error{
                .code = ErrorCode::OutOfRange,
                .message = "more than " + std::to_string(kMaxCount) + " " +
                           std::string{symbol.value()} + " atoms",
            };
        }
    }
    return formula;
}

int Formula::count(std::string_view element) const noexcept {
    const auto it = counts_.find(element);
    return it == counts_.end() ? 0 : it->second;
}

std::string Formula::hill() const {
    std::string out;
    const bool has_carbon = count("C") > 0;
    if (has_carbon) {
        append_term(out, "C", count("C"));
        if (count("H") > 0) {
            append_term(out, "H", count("H"));
        }
    }
    for (const auto& [symbol, n] : counts_) {  // std::map: alphabetical
        if (has_carbon && (symbol == "C" || symbol == "H")) {
            continue;
        }
        append_term(out, symbol, n);
    }
    return out;
}

Result<int> Formula::degrees_of_unsaturation() const {
    int twice = 2;  // 2 * DoU = 2 + sum(n * (valence - 2))
    for (const auto& [symbol, n] : counts_) {
        const auto* element = find_element(symbol);
        twice += n * (element->valence - 2);
    }
    if (twice < 0 || twice % 2 != 0) {
        return Error{
            .code = ErrorCode::InvalidArgument,
            .message = hill() +
                       " is not a closed-shell neutral molecule: degrees of "
                       "unsaturation would be " +
                       std::to_string(twice) + "/2",
        };
    }
    return twice / 2;
}

}  // namespace arise::core
