// hsqc::classify_protons and suggest_exchangeable_count on nmrshiftdb2 records with an
// assigned OH/NH/SH/SeH 1H peak (tests/data/reference/nmrshiftdb2_exchangeable.tsv, built
// by scripts/build_reference_sets.py). SECONDARY evidence: only the 1H shifts and which
// atom each H is on are measured; the HSQC is derived from the assignments (one cross peak
// per 1H peak on carbon) and the D2O-shake list here is the record's peaks on carbon. So
// these tests do not show that C-H peaks find a partner (BMRB does, bmrb_test.cpp). They
// show, on real shift distributions, how often a real exchangeable peak is hidden within
// the tolerance of another peak, that the D2O path finds every unhidden one, and that the
// suggested count never exceeds the structure's exchangeable hydrogens.

#include "arise/chem/molecule.hpp"
#include "arise/chem/perception.hpp"
#include "arise/hsqc/exchangeable.hpp"
#include "arise/io/number.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using arise::core::Ppm;
using arise::hsqc::ProtonVerdict;

// The share of exchangeable 1H peaks (mean over records) this test measured as found
// (NoHsqcPartner or Undetermined) without a D2O list at the default tolerance: 0.909722
// over 150 records on 2026-11-07, with a standard deviation of the per-record shares of
// 0.256006. Recorded rounded down and up, so the pass mark (rate less three standard
// errors) only loosens. About 9% of real exchangeable peaks lie within 0.15 ppm of a C-H or
// solvent peak and are hidden: the trade-off of the chemist's tolerance. Re-measure on
// purpose when the data or the tolerance change.
constexpr double kMeasuredFoundRate = 0.909;
constexpr double kMeasuredRecordSd = 0.257;

std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    while (true) {
        const auto end = text.find(separator);
        parts.push_back(text.substr(0, end));
        if (end == std::string_view::npos) {
            return parts;
        }
        text = text.substr(end + 1);
    }
}

double number(std::string_view text) {
    const auto value = arise::io::parse_decimal(text);
    EXPECT_TRUE(value) << text;
    return value ? value.value() : 0.0;
}

struct Peak {
    Ppm shift;
    bool exchangeable = false;  // kind X: H on O, N, S or Se
};

struct Record {
    std::string id;
    arise::nmr1d::Solvent solvent = arise::nmr1d::Solvent::CDCl3;
    std::string smiles;
    std::vector<Peak> peaks;
    std::vector<arise::io::CorrelationPeak> hsqc;
};

std::optional<arise::nmr1d::Solvent> solvent_named(std::string_view name) {
    for (const auto solvent : arise::nmr1d::kAllSolvents) {
        if (arise::nmr1d::to_string(solvent) == name) {
            return solvent;
        }
    }
    return std::nullopt;
}

Record read_record(std::string_view line) {
    const auto fields = split(line, '\t');
    EXPECT_EQ(fields.size(), 5U) << line;
    if (fields.size() != 5) {
        return {};
    }
    const auto solvent = solvent_named(fields[1]);
    EXPECT_TRUE(solvent) << fields[1];
    Record record{
        .id = std::string{fields[0]},
        .solvent = solvent.value_or(arise::nmr1d::Solvent::CDCl3),
        .smiles = std::string{fields[2]},
        .peaks = {},
        .hsqc = {},
    };
    for (const auto item : split(fields[3], ' ')) {
        const auto colon = item.find(':');
        record.peaks.push_back({
            .shift = Ppm{number(item.substr(0, colon))},
            .exchangeable = item.substr(colon + 1) == "X",
        });
    }
    if (fields[4] != "-") {
        for (const auto item : split(fields[4], ' ')) {
            const auto colon = item.find(':');
            record.hsqc.push_back({
                .proton = Ppm{number(item.substr(0, colon))},
                .carbon = Ppm{number(item.substr(colon + 1))},
                .intensity = std::nullopt,
            });
        }
    }
    return record;
}

const std::vector<Record>& records() {
    static const std::vector<Record> all = [] {
        std::vector<Record> found;
        const std::string text =
            arise::test::join_chunks(arise::test::kNmrshiftdbExchangeableChunks);
        for (const auto line : split(text, '\n')) {
            if (!line.empty() && line.front() != '#') {
                found.push_back(read_record(line));
            }
        }
        return found;
    }();
    return all;
}

std::vector<arise::io::ObservedPeak> protons(const std::vector<Peak>& peaks, bool carbon_only) {
    std::vector<arise::io::ObservedPeak> list;
    for (const auto& peak : peaks) {
        if (!carbon_only || !peak.exchangeable) {
            list.push_back({.shift = peak.shift, .integral = std::nullopt});
        }
    }
    return list;
}

arise::hsqc::ProtonClassification classify(const Record& record, bool with_d2o) {
    const auto all = protons(record.peaks, false);
    const auto after = protons(record.peaks, true);
    const auto result = arise::hsqc::classify_protons({
        .proton = all,
        .hsqc = record.hsqc,
        .solvent = record.solvent,
        .after_d2o = with_d2o ? std::optional{std::span<const arise::io::ObservedPeak>{after}}
                              : std::nullopt,
    });
    EXPECT_TRUE(result) << record.id << ": " << (result ? "" : result.error().message);
    return result ? result.value() : arise::hsqc::ProtonClassification{};
}

TEST(HsqcNmrshiftdb2, EveryRecordReadsAndHasAnExchangeablePeak) {
    ASSERT_GT(records().size(), 100U);
    for (const auto& record : records()) {
        EXPECT_TRUE(arise::chem::parse_smiles(record.smiles)) << record.id;
        EXPECT_TRUE(std::any_of(record.peaks.begin(), record.peaks.end(), [](const Peak& peak) {
            return peak.exchangeable;
        })) << record.id;
    }
}

// The share of the record's exchangeable peaks found without a D2O list; prints the others.
double found_share(const Record& record) {
    const auto classified = classify(record, false);
    int found = 0;
    int total = 0;
    for (std::size_t i = 0; i < record.peaks.size() && i < classified.peaks.size(); ++i) {
        if (!record.peaks[i].exchangeable) {
            continue;
        }
        ++total;
        const auto verdict = classified.peaks[i].verdict;
        if (verdict == ProtonVerdict::NoHsqcPartner || verdict == ProtonVerdict::Undetermined) {
            ++found;
        } else {
            std::cout << "  hidden: " << record.id << " " << record.peaks[i].shift.value()
                      << " ppm (" << arise::hsqc::to_string(verdict) << ")\n";
        }
    }
    return total == 0 ? 0.0 : static_cast<double>(found) / static_cast<double>(total);
}

TEST(HsqcNmrshiftdb2, ExchangeablePeaksAreFoundWithoutAD2oListAtTheMeasuredRate) {
    double sum = 0.0;
    double square_sum = 0.0;
    for (const auto& record : records()) {
        const double share = found_share(record);
        sum += share;
        square_sum += share * share;
    }
    const auto n = static_cast<double>(records().size());
    const double rate = sum / n;
    const double sd = std::sqrt(std::max(0.0, (square_sum - (n * rate * rate)) / (n - 1.0)));
    const double pass_mark = kMeasuredFoundRate - (3.0 * kMeasuredRecordSd / std::sqrt(n));
    std::cout << "  records " << records().size() << ", exchangeable peaks found (mean over "
              << "records) " << rate << ", SD over records " << sd << ", pass mark " << pass_mark
              << "\n";
    EXPECT_GE(rate, pass_mark);
}

// Whether another peak of the record lies within tolerance of peak `i` (so a D2O shake
// cannot show it gone).
bool overlapped(const Record& record, std::size_t i) {
    for (std::size_t j = 0; j < record.peaks.size(); ++j) {
        if (j != i && !record.peaks[j].exchangeable &&
            std::abs(record.peaks[j].shift.value() - record.peaks[i].shift.value()) <=
                arise::hsqc::kDefaultProtonTolerance.value()) {
            return true;
        }
    }
    return false;
}

// "" when, with the D2O-shake list, every unhidden exchangeable peak is shown gone and no
// C-H peak is; otherwise what is wrong. Counts exchangeable peaks under a C-H peak that the
// shake cannot show gone into `hidden`.
std::string d2o_difference(const Record& record, std::size_t& hidden) {
    const auto classified = classify(record, true);
    if (classified.peaks.size() != record.peaks.size()) {
        return "not classified";
    }
    std::string differences;
    for (std::size_t i = 0; i < record.peaks.size(); ++i) {
        const bool gone = classified.peaks[i].verdict == ProtonVerdict::Exchangeable;
        const std::string at = std::to_string(record.peaks[i].shift.value());
        if (!record.peaks[i].exchangeable) {
            differences += gone ? " C-H at " + at + " shown gone" : "";
        } else if (overlapped(record, i)) {
            hidden += gone ? 0U : 1U;
        } else {
            differences += gone ? "" : " X-H at " + at + " not shown gone";
        }
    }
    return differences;
}

TEST(HsqcNmrshiftdb2, TheD2oShakeFindsEveryUnhiddenExchangeablePeakAndNoOther) {
    std::size_t hidden = 0;
    for (const auto& record : records()) {
        EXPECT_EQ(d2o_difference(record, hidden), "") << record.id;
    }
    std::cout << "  exchangeable peaks under a C-H peak, not shown gone by the shake: " << hidden
              << "\n";
}

// "" when the count suggested with the D2O-shake list has the shown-gone peaks as its
// minimum, no more than the structure's exchangeable hydrogens, and an open maximum
// exactly in exchanging solvents; otherwise what is wrong.
std::string count_difference(const Record& record) {
    const auto classified = classify(record, true);
    const auto count = arise::hsqc::suggest_exchangeable_count(classified.peaks, record.solvent);
    const auto molecule = arise::chem::parse_smiles(record.smiles);
    if (!count || !molecule) {
        return "no count or structure";
    }
    const auto shown_gone =
        std::count_if(classified.peaks.begin(), classified.peaks.end(),
                      [](const auto& peak) { return peak.verdict == ProtonVerdict::Exchangeable; });
    const int structure = arise::chem::exchangeable_hydrogens(molecule.value());
    const bool exchanging = arise::nmr1d::solvent_peaks(record.solvent).value().exchanges_labile_h;
    std::string differences;
    if (count.value().min != shown_gone) {
        differences += " min " + std::to_string(count.value().min) + " but " +
                       std::to_string(shown_gone) + " shown gone";
    }
    if (count.value().min > structure) {
        differences += " min " + std::to_string(count.value().min) + " above the structure's " +
                       std::to_string(structure);
    }
    if (count.value().max.has_value() == exchanging) {
        differences += exchanging ? " max bounded in an exchanging solvent" : " max open";
    }
    return differences;
}

TEST(HsqcNmrshiftdb2, TheSuggestedCountNeverExceedsTheStructuresExchangeableHydrogens) {
    for (const auto& record : records()) {
        EXPECT_EQ(count_difference(record), "") << record.id;
    }
}

}  // namespace
