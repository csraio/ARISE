#include "csv_table.hpp"

#include "arise/io/number.hpp"

#include <string>
#include <utility>

namespace arise::io::detail {
namespace {

constexpr std::string_view kByteOrderMark = "\xEF\xBB\xBF";

// A line that holds data or the header, with its 1-based physical line number.
struct Line {
    std::size_t number;
    std::string_view text;
};

std::string quoted(std::string_view text) {
    return "'" + std::string{text} + "'";
}

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

// Splits on LF (dropping one trailing CR) and keeps only header and data lines:
// blank lines and `#` comments are skipped but still counted.
std::vector<Line> content_lines(std::string_view text) {
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

std::vector<std::string_view> split_fields(std::string_view line) {
    std::vector<std::string_view> fields;
    std::size_t start = 0;
    while (true) {
        const auto comma = line.find(',', start);
        fields.push_back(trim(line.substr(start, comma - start)));
        if (comma == std::string_view::npos) {
            return fields;
        }
        start = comma + 1;
    }
}

std::string column_list(std::span<const ColumnSpec> columns) {
    std::string list;
    for (const auto& column : columns) {
        list += (list.empty() ? "" : ", ") + std::string{column.name};
    }
    return list;
}

// Header field index of each spec column, if present; plus the header's field count.
struct Layout {
    std::size_t field_count;
    std::vector<std::optional<std::size_t>> positions;
};

core::Result<Layout> parse_header(const Line& header, std::span<const ColumnSpec> columns) {
    const auto names = split_fields(header.text);
    Layout layout{
        .field_count = names.size(),
        .positions = std::vector<std::optional<std::size_t>>(columns.size()),
    };
    for (std::size_t field = 0; field < names.size(); ++field) {
        const auto name = names[field];
        if (name.empty()) {
            return error_at(core::ErrorCode::ParseError, header.number,
                            "empty column name in column " + std::to_string(field + 1));
        }
        std::optional<std::size_t> spec;
        for (std::size_t i = 0; i < columns.size(); ++i) {
            if (columns[i].name == name) {
                spec = i;
            }
        }
        if (!spec.has_value()) {
            return error_at(
                core::ErrorCode::ParseError, header.number,
                "unknown column " + quoted(name) + " (allowed: " + column_list(columns) + ")");
        }
        auto& position = layout.positions[*spec];
        if (position.has_value()) {
            return error_at(core::ErrorCode::ParseError, header.number,
                            "duplicate column " + quoted(name));
        }
        position = field;
    }
    for (std::size_t i = 0; i < columns.size(); ++i) {
        if (columns[i].required && !layout.positions[i].has_value()) {
            return error_at(core::ErrorCode::ParseError, header.number,
                            "missing required column " + quoted(columns[i].name));
        }
    }
    return layout;
}

core::Result<double> parse_field(const Line& line, std::string_view column,
                                 std::string_view field) {
    if (field.empty()) {
        return error_at(core::ErrorCode::ParseError, line.number,
                        "empty " + std::string{column} + " field");
    }
    if (field.find('"') != std::string_view::npos) {
        return error_at(
            core::ErrorCode::ParseError, line.number,
            std::string{column} + " field " + quoted(field) + ": quoted fields are not supported");
    }
    auto value = parse_decimal(field);
    if (!value.ok()) {
        return error_at(core::ErrorCode::ParseError, line.number,
                        std::string{column} + ": " + value.error().message);
    }
    return value;
}

core::Result<TableRow> parse_row(const Line& line, const Layout& layout,
                                 std::span<const ColumnSpec> columns) {
    const auto fields = split_fields(line.text);
    if (fields.size() != layout.field_count) {
        return error_at(core::ErrorCode::ParseError, line.number,
                        "expected " + std::to_string(layout.field_count) +
                            (layout.field_count == 1 ? " field" : " fields") + ", got " +
                            std::to_string(fields.size()));
    }
    TableRow row{
        .line = line.number,
        .values = std::vector<std::optional<double>>(columns.size()),
        .texts = std::vector<std::string_view>(columns.size()),
    };
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const auto position = layout.positions[i];
        if (!position.has_value()) {
            continue;
        }
        const auto text = fields[*position];
        auto value = parse_field(line, columns[i].name, text);
        if (!value.ok()) {
            return std::move(value).error();
        }
        row.values[i] = value.value();
        row.texts[i] = text;
    }
    return row;
}

}  // namespace

core::Error error_at(core::ErrorCode code, std::size_t line, std::string_view what) {
    return core::Error{
        .code = code,
        .message = "line " + std::to_string(line) + ": " + std::string{what},
    };
}

core::Result<std::vector<TableRow>> read_numeric_table(std::string_view text,
                                                       std::span<const ColumnSpec> columns,
                                                       EmptyTable empty) {
    if (text.starts_with(kByteOrderMark)) {
        text.remove_prefix(kByteOrderMark.size());
    }
    const auto lines = content_lines(text);
    if (lines.empty()) {
        return core::Error{
            .code = core::ErrorCode::ParseError,
            .message = "missing header line (columns: " + column_list(columns) + ")",
        };
    }

    auto layout = parse_header(lines.front(), columns);
    if (!layout.ok()) {
        return std::move(layout).error();
    }

    std::vector<TableRow> rows;
    rows.reserve(lines.size() - 1);
    for (const Line& line : std::span{lines}.subspan(1)) {
        auto row = parse_row(line, layout.value(), columns);
        if (!row.ok()) {
            return std::move(row).error();
        }
        rows.push_back(std::move(row).value());
    }
    if (rows.empty() && empty == EmptyTable::Reject) {
        return error_at(core::ErrorCode::ParseError, lines.front().number,
                        "no peaks after the header");
    }
    return rows;
}

}  // namespace arise::io::detail
