#include "arise/io/carbon_csv.hpp"

#include "csv_table.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <utility>

namespace arise::io {
namespace {

// Indices into each reader's column table.
constexpr std::size_t kShift = 0;
constexpr std::size_t kIntensity = 1;

constexpr std::array kCarbonColumns{
    detail::ColumnSpec{.name = "shift_ppm", .required = true},
    detail::ColumnSpec{.name = "intensity", .required = false},
};
constexpr std::array kDeptColumns{
    detail::ColumnSpec{.name = "shift_ppm", .required = true},
    detail::ColumnSpec{.name = "intensity", .required = true},
};

constexpr detail::ShiftWindow kCarbonWindow{
    .nucleus = "13C",
    .min = kMinCarbonShift,
    .max = kMaxCarbonShift,
    .hint = "is it a value in Hz?",
};

core::Result<core::Ppm> carbon_shift(const detail::TableRow& row) {
    return detail::shift_in_window(row, kShift, "shift_ppm", kCarbonWindow);
}

// Converts every row with `convert`, stopping at the first error.
template <typename Peak, typename Convert>
core::Result<std::vector<Peak>> convert_rows(core::Result<std::vector<detail::TableRow>> rows,
                                             Convert convert) {
    if (!rows.ok()) {
        return std::move(rows).error();
    }
    std::vector<Peak> peaks;
    peaks.reserve(rows.value().size());
    for (const auto& row : rows.value()) {
        auto peak = convert(row);
        if (!peak.ok()) {
            return std::move(peak).error();
        }
        peaks.push_back(peak.value());
    }
    return peaks;
}

core::Result<CarbonPeak> to_carbon_peak(const detail::TableRow& row) {
    auto shift = carbon_shift(row);
    if (!shift.ok()) {
        return std::move(shift).error();
    }
    return CarbonPeak{.shift = shift.value(), .intensity = row.values[kIntensity]};
}

core::Result<DeptPeak> to_dept_peak(const detail::TableRow& row) {
    auto shift = carbon_shift(row);
    if (!shift.ok()) {
        return std::move(shift).error();
    }
    const auto intensity = row.values[kIntensity];
    assert(intensity.has_value());  // required column
    if (*intensity == 0.0) {
        return detail::error_at(core::ErrorCode::OutOfRange, row.line,
                                "intensity must be nonzero (its sign is the DEPT phase)");
    }
    return DeptPeak{.shift = shift.value(), .intensity = *intensity};
}

}  // namespace

core::Result<std::vector<CarbonPeak>> parse_carbon_csv(std::string_view text) {
    return convert_rows<CarbonPeak>(detail::read_numeric_table(text, kCarbonColumns),
                                    to_carbon_peak);
}

core::Result<std::vector<DeptPeak>> parse_dept_csv(std::string_view text) {
    return convert_rows<DeptPeak>(
        detail::read_numeric_table(text, kDeptColumns, detail::EmptyTable::Allow), to_dept_peak);
}

}  // namespace arise::io
