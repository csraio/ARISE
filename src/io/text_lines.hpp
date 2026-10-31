#pragma once
// Internal to arise_io: line handling shared by every text reader.

#include <cstddef>
#include <string_view>
#include <vector>

namespace arise::io::detail {

/// A line holding content, with its 1-based physical line number.
struct Line {
    std::size_t number;
    std::string_view text;
};

/// `text` without leading and trailing spaces and tabs.
[[nodiscard]] std::string_view trim(std::string_view text) noexcept;

/// Splits UTF-8 text into content lines: skips a leading byte-order mark, splits on LF
/// (dropping one trailing CR), and skips blank lines and lines whose first non-blank
/// character is `#`. Skipped lines still count, so numbers match the user's editor.
/// The returned views point into `text`.
[[nodiscard]] std::vector<Line> content_lines(std::string_view text);

}  // namespace arise::io::detail
