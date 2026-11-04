// The BMRB reference cases (tests/data/reference/bmrb_cases.tsv and bmrb_peaks.tsv, built
// by scripts/build_bmrb_reference.py; see BMRB_NOTICE.md) against nmr1d's solvents, and the
// notice against the cases. The cases are chosen by the builder's rule, not by hand.

#include "arise/nmr1d/solvent.hpp"
#include "bmrb_cases.hpp"
#include "reference_fixtures.hpp"

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <string_view>

namespace {

using arise::test::bmrb_cases;
namespace detail = arise::test::bmrb_detail;

TEST(BmrbReference, EveryCaseIsInAnAriseSolvent) {
    std::set<std::string_view> names;
    for (const auto solvent : arise::nmr1d::kAllSolvents) {
        names.insert(arise::nmr1d::to_string(solvent));
    }
    for (const auto& bmrb : bmrb_cases()) {
        EXPECT_TRUE(names.contains(bmrb.solvent)) << bmrb.id() << " " << bmrb.solvent;
    }
}

TEST(BmrbReference, EveryCaseHasItsOneDimensionalProtonAndHsqcLists) {
    for (const auto& bmrb : bmrb_cases()) {
        EXPECT_FALSE(bmrb.list("1H").empty()) << bmrb.id();
        EXPECT_FALSE(bmrb.list("HSQC").empty()) << bmrb.id();
    }
}

TEST(BmrbReference, NoticeListsEveryEntryOnce) {
    // BMRB asks that its data keep their authors attached: BMRB_NOTICE.md must name every
    // entry the cases use, and no other.
    const std::string notice = arise::test::join_chunks(arise::test::kBmrbNoticeChunks);
    std::multiset<std::string> listed;
    for (const auto line : detail::split(notice, '\n')) {
        if (line.starts_with("| bmse")) {
            listed.emplace(line.substr(2, line.find(' ', 2) - 2));
        }
    }
    std::set<std::string> used;
    for (const auto& bmrb : bmrb_cases()) {
        used.insert(bmrb.entry);
    }
    EXPECT_EQ(listed, std::multiset<std::string>(used.begin(), used.end()));
    for (const auto* statement : {
             "worldwide in the public domain",
             "publicly available free of charge",
             "requested to acknowledge BMRB",
             "doi: 10.1093/nar/gkac1050",
         }) {
        EXPECT_NE(notice.find(statement), std::string::npos) << statement;
    }
}

}  // namespace
