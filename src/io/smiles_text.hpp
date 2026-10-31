#pragma once
// Internal to arise_io: the character-level SMILES check shared by the readers.

#include <optional>
#include <string>
#include <string_view>

namespace arise::io::detail {

/// The first character of `smiles` outside the SMILES alphabet (ASCII letters, digits
/// and ( ) [ ] = # $ : / \ . + - @ % *), or nothing when every character is allowed.
[[nodiscard]] std::optional<char> first_non_smiles_char(std::string_view smiles) noexcept;

/// How to name a rejected character in a message: "';'", "a non-ASCII byte" or
/// "control character 0x01". Never echoes raw bytes outside printable ASCII.
[[nodiscard]] std::string describe_char(char c);

}  // namespace arise::io::detail
