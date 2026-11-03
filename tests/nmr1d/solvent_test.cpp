// Residual-solvent and water 1H peaks. Engineering cases use small hand-made tables that
// make no chemical claim; the shipped table (data/solvent_peaks.tsv) is held to every value
// the three published tables give (tests/data/reference/solvent_sources.tsv, read from the
// PDFs by scripts/build_solvent_peaks.py): each accepted shift must be confirmed by two of
// them within 0.01 ppm.
#include "arise/nmr1d/solvent.hpp"

#include "arise/core/error.hpp"
#include "arise/core/units.hpp"
#include "arise/io/number.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::core::ErrorCode;
using arise::core::Ppm;
using arise::nmr1d::kAllSolvents;
using arise::nmr1d::parse_solvent_peaks;
using arise::nmr1d::Solvent;
using arise::nmr1d::SolventPeaks;

void expect_error(std::string_view text, std::vector<std::string> fragments) {
    const auto result = parse_solvent_peaks(text);
    ASSERT_FALSE(result.ok()) << "accepted:\n" << text;
    EXPECT_EQ(result.error().code, ErrorCode::ParseError) << result.error().message;
    fragments.emplace_back("solvent peaks");
    for (const auto& fragment : fragments) {
        EXPECT_NE(result.error().message.find(fragment), std::string::npos)
            << "message '" << result.error().message << "' lacks '" << fragment << "'";
    }
}

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

// ---- reading a table ---------------------------------------------------------------------

TEST(SolventPeaks, ReadsRowsInFileOrder) {
    const auto result =
        parse_solvent_peaks("# comment\n\nCDCl3\t11\t7.2\t1.5\tno\r\nD2O\t22\t4.7,4.8\t-\tyes\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    const std::vector<SolventPeaks> expected{
        {
            .solvent = Solvent::CDCl3,
            .pubchem_cid = 11,
            .residual = {Ppm{7.2}},
            .water = {Ppm{1.5}},
            .exchanges_labile_h = false,
        },
        {
            .solvent = Solvent::D2O,
            .pubchem_cid = 22,
            .residual = {Ppm{4.7}, Ppm{4.8}},
            .water = {},
            .exchanges_labile_h = true,
        },
    };
    EXPECT_EQ(result.value(), expected);
}

TEST(SolventPeaks, AcceptsEqualShiftsAndTheWindowEdges) {
    const auto result = parse_solvent_peaks("THF-d8\t1\t-10,25,25\t-10\tno\n");
    ASSERT_TRUE(result.ok()) << result.error().message;
    EXPECT_EQ(result.value()[0].residual, (std::vector{Ppm{-10}, Ppm{25}, Ppm{25}}));
}

TEST(SolventPeaks, RejectsMalformedRows) {
    expect_error("", {"no rows"});
    expect_error("# only comments\n", {"no rows"});
    expect_error("CDCl3\t1\t7.2\t1.5\n", {"line 1", "expected 5 fields"});
    expect_error("CDCl3\t1\t7.2\t1.5\tno\textra\n", {"line 1", "expected 5 fields"});
    expect_error("# c\nchloroform\t1\t7.2\t1.5\tno\n",
                 {"line 2", "unknown solvent 'chloroform'", "CDCl3", "(CD3)2SO", "D2O"});
    expect_error("cdcl3\t1\t7.2\t1.5\tno\n", {"unknown solvent 'cdcl3'"});
    expect_error("CDCl3\t1\t7.2\t1.5\tno\nCDCl3\t2\t7.3\t1.6\tno\n",
                 {"line 2", "CDCl3 is listed twice"});
    expect_error("CDCl3\t0\t7.2\t1.5\tno\n", {"line 1", "CID", "'0'"});
    expect_error("CDCl3\tx\t7.2\t1.5\tno\n", {"line 1", "CID", "'x'"});
    expect_error("CDCl3\t-1\t7.2\t1.5\tno\n", {"line 1", "CID", "'-1'"});
    expect_error("CDCl3\t1\tyes\t1.5\tno\n", {"line 1", "not a number: 'yes'"});
    expect_error("CDCl3\t1\t-\t1.5\tno\n", {"line 1", "no residual shift"});
    expect_error("CDCl3\t1\t7.2,\t1.5\tno\n", {"line 1", "not a number: ''"});
    expect_error("CDCl3\t1\t7.2\t1.5,,1.6\tno\n", {"line 1", "not a number: ''"});
    expect_error("CDCl3\t1\t7.2\t\tno\n", {"line 1", "not a number: ''"});
    expect_error("CDCl3\t1\t7.2\t1.5\tYes\n", {"line 1", "'yes' or 'no'", "'Yes'"});
    expect_error("CDCl3\t1\t7.2\t1.5\t\n", {"line 1", "'yes' or 'no'"});
}

TEST(SolventPeaks, RejectsShiftsOutsideTheProtonWindowOrOutOfOrder) {
    expect_error("CDCl3\t1\t25.01\t1.5\tno\n", {"line 1", "25.01", "outside"});
    expect_error("CDCl3\t1\t7.2\t-10.01\tno\n", {"line 1", "-10.01", "outside"});
    expect_error("CDCl3\t1\t7.2\tinf\tno\n", {"line 1", "not a number: 'inf'"});
    expect_error("CDCl3\t1\t7.3,7.2\t1.5\tno\n", {"line 1", "7.2", "ascending"});
    expect_error("CDCl3\t1\t7.2\t1.6,1.5\tno\n", {"line 1", "1.5", "ascending"});
}

TEST(SolventPeaks, NamesAreDistinct) {
    std::set<std::string_view> names;
    for (const auto solvent : kAllSolvents) {
        EXPECT_FALSE(to_string(solvent).empty());
        EXPECT_TRUE(names.insert(to_string(solvent)).second) << to_string(solvent);
    }
}

// ---- the shipped table -------------------------------------------------------------------

TEST(SolventPeaks, BuiltInTableListsEverySolventOnceUnderItsPrintedName) {
    // parse_solvent_peaks matches names exactly, so this also holds every to_string to the
    // heading the generated table prints.
    const auto table =
        parse_solvent_peaks(arise::test::join_chunks(arise::test::kSolventPeaksChunks));
    ASSERT_TRUE(table.ok()) << table.error().message;
    std::multiset<Solvent> listed;
    for (const auto& row : table.value()) {
        listed.insert(row.solvent);
    }
    EXPECT_EQ(listed, std::multiset<Solvent>(kAllSolvents.begin(), kAllSolvents.end()));
    for (const auto& row : table.value()) {
        const auto built_in = arise::nmr1d::solvent_peaks(row.solvent);
        ASSERT_TRUE(built_in.ok()) << built_in.error().message;
        EXPECT_EQ(built_in.value(), row);
    }
}

// One row of solvent_sources.tsv: the shifts a source (or "accepted", or "one source only")
// gives for one kind of peak of one solvent.
struct SourceRow {
    std::string solvent;
    std::string kind;
    std::string source;
    std::vector<double> shifts;
};

using Key = std::pair<std::string, std::string>;  // solvent, kind

bool is_publication(const SourceRow& row) {
    return row.source != "accepted" && row.source != "one source only";
}

std::vector<double> shift_list(std::string_view field) {
    std::vector<double> shifts;
    if (field == "-") {
        return shifts;
    }
    for (const auto item : split(field, ',')) {
        const auto value = arise::io::parse_decimal(item);
        EXPECT_TRUE(value.ok()) << field;
        shifts.push_back(value.value_or(0.0));
    }
    return shifts;
}

std::vector<SourceRow> source_rows() {
    std::vector<SourceRow> rows;
    // Kept alive for the loop: split returns views into it.
    const auto text = arise::test::join_chunks(arise::test::kSolventSourcesChunks);
    for (const auto line : split(text, '\n')) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const auto fields = split(line, '\t');
        EXPECT_EQ(fields.size(), 6U) << line;
        if (fields.size() == 6U) {
            rows.push_back({
                .solvent = std::string{fields[0]},
                .kind = std::string{fields[2]},
                .source = std::string{fields[3]},
                .shifts = shift_list(fields[5]),
            });
        }
    }
    return rows;
}

// The rows of one source label ("accepted" or "one source only"), by solvent and kind.
std::map<Key, std::vector<double>> rows_labelled(const std::vector<SourceRow>& rows,
                                                 std::string_view label) {
    std::map<Key, std::vector<double>> found;
    for (const auto& row : rows) {
        if (row.source == label) {
            found[{row.solvent, row.kind}] = row.shifts;
        }
    }
    return found;
}

constexpr double kAgreement = 0.01 + 1e-9;  // the builder's rule; slack for binary rounding

bool within(double a, double b) {
    return std::fabs(a - b) <= kAgreement;
}

std::vector<double> values_of(const std::vector<Ppm>& shifts) {
    std::vector<double> values;
    values.reserve(shifts.size());
    for (const auto shift : shifts) {
        values.push_back(shift.value());
    }
    return values;
}

// The publications giving a value within 0.01 ppm of `shift` for this solvent and kind.
std::set<std::string> confirming(const std::vector<SourceRow>& rows, const Key& key, double shift) {
    std::set<std::string> sources;
    for (const auto& row : rows) {
        const bool counts = is_publication(row) && Key{row.solvent, row.kind} == key;
        for (const double value : counts ? row.shifts : std::vector<double>{}) {
            if (within(value, shift)) {
                sources.insert(row.source);
            }
        }
    }
    return sources;
}

// The shipped residual and water shifts of every solvent, by solvent and kind.
std::map<Key, std::vector<double>> shipped_shifts() {
    std::map<Key, std::vector<double>> shipped;
    for (const auto solvent : kAllSolvents) {
        const auto peaks = arise::nmr1d::solvent_peaks(solvent);
        EXPECT_TRUE(peaks.ok()) << to_string(solvent);
        const std::string name{to_string(solvent)};
        shipped[{name, "residual"}] = values_of(peaks.value_or(SolventPeaks{}).residual);
        shipped[{name, "water"}] = values_of(peaks.value_or(SolventPeaks{}).water);
    }
    return shipped;
}

TEST(SolventPeaks, EveryShippedShiftIsConfirmedByTwoPublishedTables) {
    const auto rows = source_rows();
    ASSERT_FALSE(rows.empty());
    for (const auto& [key, shifts] : shipped_shifts()) {
        for (const double shift : shifts) {
            EXPECT_GE(confirming(rows, key, shift).size(), 2U)
                << key.first << " " << key.second << " " << shift;
        }
    }
}

TEST(SolventPeaks, EverySourceValueIsShippedOrListedAsUnconfirmed) {
    // Nothing a table gives is silently dropped: it is within 0.01 ppm of a shipped shift,
    // or it is listed under "one source only" for the chemist to see.
    const auto rows = source_rows();
    auto accepted = rows_labelled(rows, "accepted");
    auto unconfirmed = rows_labelled(rows, "one source only");
    for (const auto& row : rows) {
        const Key key{row.solvent, row.kind};
        for (const double value : is_publication(row) ? row.shifts : std::vector<double>{}) {
            const auto& near = accepted[key];
            const auto& alone = unconfirmed[key];
            const bool kept =
                std::any_of(near.begin(), near.end(), [&](double s) { return within(value, s); }) ||
                std::find(alone.begin(), alone.end(), value) != alone.end();
            EXPECT_TRUE(kept) << row.solvent << " " << row.kind << " " << row.source << " "
                              << value;
        }
    }
}

TEST(SolventPeaks, ShippedTableIsExactlyTheAcceptedValues) {
    // Every solvent with an accepted residual peak ships, with exactly the accepted shifts.
    const auto accepted = rows_labelled(source_rows(), "accepted");
    std::map<Key, std::vector<double>> expected;
    for (const auto& [key, shifts] : accepted) {
        const auto residual = accepted.find({key.first, "residual"});
        if (residual != accepted.end() && !residual->second.empty()) {
            expected[key] = shifts;
        }
    }
    EXPECT_EQ(shipped_shifts(), expected);
}

}  // namespace
