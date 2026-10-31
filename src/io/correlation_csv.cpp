#include "arise/io/correlation_csv.hpp"

#include "arise/io/carbon_csv.hpp"
#include "arise/io/peak_csv.hpp"
#include "csv_table.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace arise::io {
namespace {

// Indices into kColumns.
constexpr std::size_t kProton = 0;
constexpr std::size_t kCarbon = 1;
constexpr std::size_t kIntensity = 2;

constexpr std::array kColumns{
    detail::ColumnSpec{.name = "h_ppm", .required = true},
    detail::ColumnSpec{.name = "c_ppm", .required = true},
    detail::ColumnSpec{.name = "intensity", .required = false},
};

constexpr std::string_view kSwapHint = "are the h_ppm and c_ppm columns swapped?";
constexpr detail::ShiftWindow kProtonWindow{
    .nucleus = "1H",
    .min = kMinProtonShift,
    .max = kMaxProtonShift,
    .hint = kSwapHint,
};
constexpr detail::ShiftWindow kCarbonWindow{
    .nucleus = "13C",
    .min = kMinCarbonShift,
    .max = kMaxCarbonShift,
    .hint = kSwapHint,
};

// Shared by both experiments; HSQC additionally needs a nonzero intensity.
core::Result<CorrelationPeak> to_correlation(const detail::TableRow& row) {
    auto proton = detail::shift_in_window(row, kProton, "h_ppm", kProtonWindow);
    if (!proton.ok()) {
        return std::move(proton).error();
    }
    auto carbon = detail::shift_in_window(row, kCarbon, "c_ppm", kCarbonWindow);
    if (!carbon.ok()) {
        return std::move(carbon).error();
    }
    return CorrelationPeak{
        .proton = proton.value(),
        .carbon = carbon.value(),
        .intensity = row.values[kIntensity],
    };
}

core::Result<CorrelationPeak> to_hsqc_peak(const detail::TableRow& row) {
    auto peak = to_correlation(row);
    if (peak.ok() && peak.value().intensity == 0.0) {
        return detail::error_at(core::ErrorCode::OutOfRange, row.line,
                                "intensity must be nonzero (its sign is the HSQC edit phase)");
    }
    return peak;
}

}  // namespace

core::Result<std::vector<CorrelationPeak>> parse_hsqc_csv(std::string_view text) {
    return detail::convert_rows<CorrelationPeak>(detail::read_numeric_table(text, kColumns),
                                                 to_hsqc_peak);
}

core::Result<std::vector<CorrelationPeak>> parse_hmbc_csv(std::string_view text) {
    return detail::convert_rows<CorrelationPeak>(detail::read_numeric_table(text, kColumns),
                                                 to_correlation);
}

}  // namespace arise::io
