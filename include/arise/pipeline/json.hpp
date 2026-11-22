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
#include "arise/nmr1d/dept.hpp"
#include "arise/pipeline/pipeline.hpp"

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

/// The classified 13C peaks, in 13C order:
/// [{"shift":29.5,"types":["CH2"]}, ...], types in C, CH, CH2, CH3 order.
[[nodiscard]] std::string to_json(std::span<const nmr1d::CarbonAssignment> peaks);

/// Group suggestions, in chem::kAllGroups order:
/// [{"group":"ketone","evidence":[{"peak":0,"shift":208.1,"hallmark":"...","types":[...],
///   "low":195.2,"high":215.7,"n":812,"molecules":790,"source":"percentiles"}]}]
/// ("source": "percentiles", "measured range" or "lab").
[[nodiscard]] std::string to_json(std::span<const filter::GroupSuggestion> suggestions);

/// {"formula":"C15H20O3","rows":13040,
///  "removed":[{"cid":7,"smiles":"...","reason":"mixture"}],
///  "candidates":[{"cid":1,"smiles":"...","tier":"active","reasons":[...],"notes":[...],
///   "assumed_overlaps":0,"environments_fewest":12,"environments_most":14}]}
/// ("tier": "active", "back-burner" or "eliminated"; "reason" as io::to_string).
[[nodiscard]] std::string to_json(const Ranking& ranking);

}  // namespace arise::pipeline
