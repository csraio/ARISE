#pragma once
// JSON text for the web client: the page reads every pipeline result with JavaScript's
// JSON.parse. Written in house (no JSON library). Object keys come in a fixed order and
// numbers in a fixed form, so equal results give equal text, natively and in WebAssembly.
//
// Shifts are written like the engine's messages: fixed-point, six decimals at most,
// trailing zeros dropped (29.5, 170). Counts are whole numbers. Strings are UTF-8 with
// '"', '\' and control characters escaped.

#include "arise/core/error.hpp"
#include "arise/filter/groups.hpp"
#include "arise/filter/proton_groups.hpp"
#include "arise/io/peak_csv.hpp"
#include "arise/nmr1d/dept.hpp"
#include "arise/pipeline/pipeline.hpp"

#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace arise::pipeline {

/// `text` as a JSON string, quotes included.
[[nodiscard]] std::string json_string(std::string_view text);

/// {"error":{"code":"ParseError","message":"..."}} (core::to_string names the code).
[[nodiscard]] std::string to_json(const core::Error& error);

/// {"url":"..."}
[[nodiscard]] std::string url_json(std::string_view url);

/// Every group's name (chem::to_string), in chem::kAllGroups order: ["carboxylic acid", ...].
[[nodiscard]] std::string groups_json();

/// The classified 13C peaks, in 13C order:
/// [{"shift":29.5,"types":["CH2"]}, ...], types in C, CH, CH2, CH3 order.
[[nodiscard]] std::string to_json(std::span<const nmr1d::CarbonAssignment> peaks);

/// Group suggestions of both spectra, merged per group in chem::kAllGroups order, 13C
/// evidence before 1H evidence:
/// [{"group":"methoxy","evidence":[
///   {"spectrum":"13C","peak":0,"shift":55.9,"hallmark":"...","types":[...],"low":49.3,
///    "high":62.7,"n":3423,"molecules":2234,"source":"percentiles"},
///   {"spectrum":"1H","peak":2,"shift":3.81,"protons":3,"hallmark":"O-CH3 carbon",
///    "h_count":"CH3","low":3.04,"high":4.31,"n":1082,"molecules":807,"precision":0.719,
///    "base_rate":0.14}]}]
/// ("source": "percentiles", "measured range" or "lab"; "protons" null without an integral;
/// rates with at most three decimals).
[[nodiscard]] std::string to_json(const Suggestions& suggestions);

/// The 1H peaks: [{"shift":3.81,"protons":3}, ...] ("protons" as filter::proton_count, null
/// without an integral).
[[nodiscard]] std::string to_json(std::span<const io::ObservedPeak> peaks);

/// to_json of the 1H peaks, or null when no 1H spectrum was entered.
[[nodiscard]] std::string protons_json(const std::optional<std::vector<io::ObservedPeak>>& peaks);

/// {"formula":"C15H20O3","rows":13040,
///  "removed":[{"cid":7,"smiles":"...","reason":"mixture"}],
///  "candidates":[{"cid":1,"smiles":"...","tier":"active","reasons":[...],"notes":[...],
///   "assumed_overlaps":0,"environments_fewest":12,"environments_most":14}]}
/// ("tier": "active", "back-burner" or "eliminated"; "reason" as io::to_string).
[[nodiscard]] std::string to_json(const Ranking& ranking);

}  // namespace arise::pipeline
