#include "arise/io/window_csv.hpp"

#include "arise/io/carbon_csv.hpp"
#include "csv_table.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <utility>

namespace arise::io {
namespace {

// Indices into kWindowColumns.
constexpr std::size_t kGroup = 0;
constexpr std::size_t kHallmark = 1;
constexpr std::size_t kDept = 2;
constexpr std::size_t kLow = 3;
constexpr std::size_t kHigh = 4;

constexpr std::array kWindowColumns{
    detail::ColumnSpec{.name = "group", .required = true, .kind = detail::ColumnKind::Text},
    detail::ColumnSpec{.name = "hallmark", .required = false, .kind = detail::ColumnKind::Text},
    detail::ColumnSpec{.name = "dept", .required = false, .kind = detail::ColumnKind::Text},
    detail::ColumnSpec{.name = "low_ppm", .required = true},
    detail::ColumnSpec{.name = "high_ppm", .required = true},
};

constexpr detail::ShiftWindow kCarbonWindow{
    .nucleus = "13C",
    .min = kMinCarbonShift,
    .max = kMaxCarbonShift,
    .hint = "is it a value in Hz?",
};

core::Result<WindowRow> to_window_row(const detail::TableRow& row) {
    auto low = detail::shift_in_window(row, kLow, "low_ppm", kCarbonWindow);
    if (!low.ok()) {
        return std::move(low).error();
    }
    auto high = detail::shift_in_window(row, kHigh, "high_ppm", kCarbonWindow);
    if (!high.ok()) {
        return std::move(high).error();
    }
    if (low.value() > high.value()) {
        return detail::error_at(core::ErrorCode::OutOfRange, row.line,
                                "low_ppm " + std::string{row.texts[kLow]} + " is above high_ppm " +
                                    std::string{row.texts[kHigh]});
    }
    return WindowRow{
        .line = row.line,
        .group = std::string{row.texts[kGroup]},
        .hallmark = std::string{row.texts[kHallmark]},
        .dept = std::string{row.texts[kDept]},
        .low = low.value(),
        .high = high.value(),
    };
}

}  // namespace

core::Result<std::vector<WindowRow>> parse_window_csv(std::string_view text) {
    return detail::convert_rows<WindowRow>(
        detail::read_numeric_table(text, kWindowColumns, detail::EmptyTable::Reject, "rows"),
        to_window_row);
}

}  // namespace arise::io
