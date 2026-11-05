// hsqc::classify_protons against BMRB's experimental spectra (tests/data/reference/
// bmrb_cases.tsv and bmrb_peaks.tsv, built by scripts/build_bmrb_reference.py; see
// BMRB_NOTICE.md). The cases are chosen by the builder's rule, and the expected verdicts
// come from BMRB's own assignments: which atom each 1H peak's H is bonded to, read from the
// entry's structure. The input is each case's own 1H and HSQC lists, unmodified, with the
// HSQC passed through ARISE's reader (io::parse_hsqc_csv) as a user's file would be.
//
// Checked:
//   * 1H peaks BMRB assigns to H on carbon must find an HSQC partner (OnCarbon, or
//     Undetermined when a solvent signal also lies near) at a measured rate. A per-peak
//     rule is not possible: a 1H list gives a multiplet's centre, an HSQC list a point
//     inside it, and BMRB's two lists do not always line up.
//   * 1H peaks BMRB assigns to H on O, N or S must not be found on carbon (NoHsqcPartner,
//     or Undetermined). Disagreements are signed-off rows of known_differences.tsv.

#include "arise/chem/molecule.hpp"
#include "arise/hsqc/exchangeable.hpp"
#include "arise/io/correlation_csv.hpp"
#include "bmrb_cases.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using arise::chem::Element;
using arise::hsqc::ProtonVerdict;
using arise::test::bmrb_cases;
using arise::test::BmrbCase;
namespace detail = arise::test::bmrb_detail;

// The share of BMRB C-H 1H peaks (mean over cases) that this test measured as finding an
// HSQC partner at the default tolerance: 0.984566 over 682 cases on 2026-11-05, recorded
// rounded down. The pass mark is this rate less three standard errors with cases as the
// units (the hold-out convention of tests/filter/holdout_test.cpp), so a change that makes
// matching worse fails. The statistic is a mean of per-case fractions, not a count of
// yes/no outcomes, so the standard error is the measured standard deviation of those
// fractions over sqrt(n) rather than the binomial sqrt(p(1-p)/n), which for values in
// [0, 1] is only an upper bound (here 0.0048 against the true 0.0030). Standard deviation
// measured on the same run: 0.0782431, recorded rounded up. Both roundings loosen the mark.
// Re-measure both on purpose when the data or the tolerance change.
constexpr double kMeasuredRate = 0.984;
constexpr double kMeasuredCaseSd = 0.079;

std::optional<arise::nmr1d::Solvent> solvent_named(std::string_view name) {
    for (const auto solvent : arise::nmr1d::kAllSolvents) {
        if (arise::nmr1d::to_string(solvent) == name) {
            return solvent;
        }
    }
    return std::nullopt;
}

// The case's HSQC list as a user's CSV file, read by ARISE's reader.
arise::core::Result<std::vector<arise::io::CorrelationPeak>> read_hsqc(const BmrbCase& bmrb) {
    std::string text = "h_ppm,c_ppm\n";
    for (const auto& peak : bmrb.list("HSQC")) {
        text += std::to_string(peak.h_shift.value_or(0.0)) + "," +
                std::to_string(peak.c_shift.value_or(0.0)) + "\n";
    }
    return arise::io::parse_hsqc_csv(text);
}

enum class Bonded : std::uint8_t { Carbon, Heteroatom, Mixed };

// What the atoms BMRB assigns a 1H peak to are, in ARISE's reading of the structure.
Bonded bonded_to(const arise::chem::Molecule& molecule, const std::vector<std::size_t>& atoms) {
    bool carbon = false;
    bool hetero = false;
    for (const auto atom : atoms) {
        const auto element = molecule.atoms()[atom].element;
        carbon = carbon || element == Element::C;
        hetero = hetero || element == Element::O || element == Element::N || element == Element::S;
    }
    if (carbon == hetero) {
        return Bonded::Mixed;
    }
    return carbon ? Bonded::Carbon : Bonded::Heteroatom;
}

struct CaseResult {
    std::size_t carbon_peaks = 0;
    std::size_t carbon_found = 0;
    std::vector<std::string> carbon_misses;
    std::vector<std::string> hetero_wrong;
};

CaseResult check_case(const BmrbCase& bmrb, const std::vector<arise::io::CorrelationPeak>& hsqc) {
    CaseResult result;
    const auto molecule = arise::chem::parse_smiles(bmrb.smiles);
    const auto solvent = solvent_named(bmrb.solvent);
    if (!molecule || !solvent) {
        result.hetero_wrong.emplace_back("structure or solvent not read");
        return result;
    }
    const auto list = bmrb.list("1H");
    std::vector<arise::io::ObservedPeak> protons;
    protons.reserve(list.size());
    for (const auto& peak : list) {
        protons.push_back(
            {.shift = arise::core::Ppm{peak.h_shift.value_or(0.0)}, .integral = std::nullopt});
    }
    const auto classified =
        arise::hsqc::classify_protons({.proton = protons, .hsqc = hsqc, .solvent = *solvent});
    if (!classified) {
        result.hetero_wrong.push_back("classify_protons: " + classified.error().message);
        return result;
    }
    for (std::size_t i = 0; i < list.size(); ++i) {
        if (list[i].h_atoms.empty()) {
            continue;
        }
        const auto verdict = classified.value().peaks[i].verdict;
        const bool partnered =
            verdict == ProtonVerdict::OnCarbon || verdict == ProtonVerdict::Undetermined;
        const std::string where = std::to_string(list[i].h_shift.value_or(0.0)) + " ppm (" +
                                  std::string{arise::hsqc::to_string(verdict)} + ")";
        switch (bonded_to(molecule.value(), list[i].h_atoms)) {
            case Bonded::Carbon:
                ++result.carbon_peaks;
                result.carbon_found += partnered ? 1U : 0U;
                if (!partnered) {
                    result.carbon_misses.push_back(where);
                }
                break;
            case Bonded::Heteroatom:
                if (verdict != ProtonVerdict::NoHsqcPartner &&
                    verdict != ProtonVerdict::Undetermined) {
                    result.hetero_wrong.push_back(where);
                }
                break;
            case Bonded::Mixed:
                break;
        }
    }
    return result;
}

// Signed-off differences of source "bmrb" in known_differences.tsv; each must still occur.
std::set<std::pair<std::string, std::string>> known_bmrb_differences() {
    std::set<std::pair<std::string, std::string>> rows;
    const std::string text = arise::test::join_chunks(arise::test::kKnownDifferencesChunks);
    for (const auto line : detail::data_lines(text)) {
        const auto fields = detail::split(line, '\t');
        if (fields.size() == 4 && fields[0] == "bmrb") {
            rows.insert({std::string{fields[1]}, std::string{fields[2]}});
        }
    }
    return rows;
}

struct Survey {
    std::size_t cases = 0;
    std::size_t rejected = 0;
    double rate_sum = 0.0;
    double rate_sq_sum = 0.0;
};

// Adds one case to the survey: its C-H rate, and a failure for each H on O/N/S found on
// carbon unless the case is a signed-off difference (removed from `unused` when it is).
void survey_case(const BmrbCase& bmrb, Survey& survey,
                 std::set<std::pair<std::string, std::string>>& unused) {
    const auto hsqc = read_hsqc(bmrb);
    if (!hsqc) {
        // BMRB lists some HSQC axes the wrong way round; ARISE's reader must refuse them as
        // out of range, as it would a user's swapped columns.
        EXPECT_EQ(hsqc.error().code, arise::core::ErrorCode::OutOfRange) << bmrb.id();
        ++survey.rejected;
        return;
    }
    const auto result = check_case(bmrb, hsqc.value());
    for (const auto& miss : result.carbon_misses) {
        std::cout << "  C-H without HSQC partner: " << bmrb.id() << " " << miss << "\n";
    }
    if (!result.hetero_wrong.empty() && unused.erase({bmrb.id(), "exchangeable"}) == 0) {
        for (const auto& wrong : result.hetero_wrong) {
            ADD_FAILURE() << bmrb.id() << " " << bmrb.compound << ": H on O/N/S found on carbon at "
                          << wrong;
        }
    }
    if (result.carbon_peaks > 0) {
        ++survey.cases;
        const double rate =
            static_cast<double>(result.carbon_found) / static_cast<double>(result.carbon_peaks);
        survey.rate_sum += rate;
        survey.rate_sq_sum += rate * rate;
    }
}

TEST(HsqcBmrb, ProtonsOnCarbonFindTheirHsqcPartnerAtTheMeasuredRate) {
    Survey survey;
    auto unused = known_bmrb_differences();
    for (const auto& bmrb : bmrb_cases()) {
        survey_case(bmrb, survey, unused);
    }
    for (const auto& [molecule, check] : unused) {
        ADD_FAILURE() << molecule << " " << check << " is listed but now agrees";
    }
    ASSERT_GT(survey.cases, 600U);
    const auto n = static_cast<double>(survey.cases);
    const double rate = survey.rate_sum / n;
    // Sample standard deviation of the per-case rates (n - 1 denominator).
    const double case_sd =
        std::sqrt(std::max(0.0, (survey.rate_sq_sum - (n * rate * rate)) / (n - 1.0)));
    const double pass_mark = kMeasuredRate - (3.0 * kMeasuredCaseSd / std::sqrt(n));
    std::cout << "  cases " << survey.cases << ", HSQC lists refused " << survey.rejected
              << ", C-H peaks with an HSQC partner (mean over cases) " << rate << ", SD over cases "
              << case_sd << ", pass mark " << pass_mark << "\n";
    EXPECT_GE(rate, pass_mark);
}

}  // namespace
