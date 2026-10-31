#include "smiles_text.hpp"

namespace arise::io::detail {
namespace {

bool is_smiles_char(char c) noexcept {
    if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')) {
        return true;
    }
    constexpr std::string_view kSymbols = "()[]=#$:/\\.+-@%*";
    return kSymbols.find(c) != std::string_view::npos;
}

}  // namespace

std::optional<char> first_non_smiles_char(std::string_view smiles) noexcept {
    for (const char c : smiles) {
        if (!is_smiles_char(c)) {
            return c;
        }
    }
    return std::nullopt;
}

std::string describe_char(char c) {
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

}  // namespace arise::io::detail
