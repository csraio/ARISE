#include "arise/io/peak_csv.hpp"

#include "csv_table.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <utility>

namespace arise::io {
namespace {

// Indices into kColumns.
constexpr std::size_t kShift = 0;
constexpr std::size_t kIntegral = 1;
constexpr std::array kColumns{
    detail::ColumnSpec{.name = "shift_ppm", .required = true},
    detail::ColumnSpec{.name = "integral", .required = false},
};

constexpr detail::ShiftWindow kProtonWindow{
    .nucleus = "1H",
    .min = kMinProtonShift,
    .max = kMaxProtonShift,
    .hint = "is it a 13C shift or a value in Hz?",
};

core::Result<ObservedPeak> to_peak(const detail::TableRow& row) {
    auto shift = detail::shift_in_window(row, kShift, "shift_ppm", kProtonWindow);
    if (!shift.ok()) {
        return std::move(shift).error();
    }
    const auto integral = row.values[kIntegral];
    if (integral.has_value() && *integral <= 0.0) {
        return detail::error_at(core::ErrorCode::OutOfRange, row.line,
                                "integral must be > 0, got " + std::string{row.texts[kIntegral]});
    }
    return ObservedPeak{.shift = shift.value(), .integral = integral};
}

}  // namespace

core::Result<std::vector<ObservedPeak>> parse_peak_csv(std::string_view text) {
    auto rows = detail::read_numeric_table(text, kColumns);
    if (!rows.ok()) {
        return std::move(rows).error();
    }
    std::vector<ObservedPeak> peaks;
    peaks.reserve(rows.value().size());
    for (const auto& row : rows.value()) {
        auto peak = to_peak(row);
        if (!peak.ok()) {
            return std::move(peak).error();
        }
        peaks.push_back(peak.value());
    }
    return peaks;
}

}  // namespace arise::io
