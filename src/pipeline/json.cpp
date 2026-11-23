#include "arise/pipeline/json.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>

namespace arise::pipeline {
namespace {

// A shift as the engine's messages write it: fixed-point, six decimals at most, trailing
// zeros dropped; "-0" becomes "0".
std::string number(core::Ppm shift) {
    std::string text = std::to_string(shift.value());
    text.erase(text.find_last_not_of('0') + 1);
    if (text.back() == '.') {
        text.pop_back();
    }
    return text == "-0" ? "0" : text;
}

std::string number(std::size_t value) {
    return std::to_string(value);
}

std::string number(int value) {
    return std::to_string(value);
}

// A JSON array of `items`, each written by `write`.
template <typename Items, typename Write>
std::string array(const Items& items, Write write) {
    std::string text = "[";
    bool first = true;
    for (const auto& item : items) {
        text += first ? "" : ",";
        first = false;
        text += write(item);
    }
    return text + "]";
}

std::string types(nmr1d::CarbonTypeSet set) {
    constexpr std::array kTypes{
        nmr1d::CarbonType::C,
        nmr1d::CarbonType::CH,
        nmr1d::CarbonType::CH2,
        nmr1d::CarbonType::CH3,
    };
    std::string text = "[";
    for (const auto type : kTypes) {
        if (set.contains(type)) {
            text += (text.size() > 1 ? "," : "") + json_string(nmr1d::to_string(type));
        }
    }
    return text + "]";
}

std::string_view source_name(filter::WindowSource source) {
    switch (source) {
        case filter::WindowSource::Percentiles:
            return "percentiles";
        case filter::WindowSource::MeasuredRange:
            return "measured range";
        case filter::WindowSource::Lab:
            return "lab";
    }
    return "";
}

// A rate with at most three decimals, trailing zeros dropped (0.719, 0.14, 1).
std::string rate(double value) {
    std::string text = std::to_string(std::round(value * 1000.0) / 1000.0);
    text.erase(text.find_last_not_of('0') + 1);
    if (text.back() == '.') {
        text.pop_back();
    }
    return text == "-0" ? "0" : text;
}

std::string optional_count(std::optional<int> value) {
    return value.has_value() ? std::to_string(*value) : "null";
}

std::string evidence(const filter::PeakEvidence& e) {
    return R"({"spectrum":"13C","peak":)" + number(e.peak) + ",\"shift\":" + number(e.shift) +
           ",\"hallmark\":" + json_string(e.window.hallmark) +
           ",\"types\":" + types(e.window.types) + ",\"low\":" + number(e.window.low) +
           ",\"high\":" + number(e.window.high) + ",\"n\":" + number(e.window.n) +
           ",\"molecules\":" + number(e.window.molecules) +
           ",\"source\":" + json_string(source_name(e.window.source)) + "}";
}

std::string proton_evidence(const filter::ProtonEvidence& e) {
    return R"({"spectrum":"1H","peak":)" + number(e.peak) + ",\"shift\":" + number(e.shift) +
           ",\"protons\":" + optional_count(e.protons) +
           ",\"hallmark\":" + json_string(e.window.hallmark) +
           ",\"h_count\":" + json_string(nmr1d::to_string(e.window.type)) +
           ",\"low\":" + number(e.window.low) + ",\"high\":" + number(e.window.high) +
           ",\"n\":" + number(e.window.n) + ",\"molecules\":" + number(e.window.molecules) +
           ",\"precision\":" + rate(e.window.precision) +
           ",\"base_rate\":" + rate(e.window.base_rate) + "}";
}

std::string strings(const std::vector<std::string>& items) {
    return array(items, [](const std::string& item) { return json_string(item); });
}

std::string candidate(const RankedCandidate& c) {
    return "{\"cid\":" + std::to_string(c.cid) + ",\"smiles\":" + json_string(c.smiles) +
           ",\"tier\":" + json_string(filter::to_string(c.verdict.tier)) +
           ",\"reasons\":" + strings(c.verdict.reasons) + ",\"notes\":" + strings(c.notes) +
           ",\"assumed_overlaps\":" + number(c.assumed_overlaps) +
           ",\"environments_fewest\":" + number(c.environments_fewest) +
           ",\"environments_most\":" + number(c.environments_most) + "}";
}

// The evidence of `group` from both spectra, 13C first, comma-separated ("" if none).
std::string group_evidence(const Suggestions& suggestions, chem::Group group) {
    std::string items;
    const auto add = [&items](const std::string& item) {
        items += (items.empty() ? "" : ",") + item;
    };
    for (const auto& s : suggestions.carbon) {
        if (s.group == group) {
            std::for_each(s.evidence.begin(), s.evidence.end(),
                          [&](const auto& e) { add(evidence(e)); });
        }
    }
    for (const auto& s : suggestions.proton) {
        if (s.group == group) {
            std::for_each(s.evidence.begin(), s.evidence.end(),
                          [&](const auto& e) { add(proton_evidence(e)); });
        }
    }
    return items;
}

}  // namespace

std::string json_string(std::string_view text) {
    constexpr std::string_view kHex = "0123456789abcdef";
    std::string quoted = "\"";
    for (const char c : text) {
        switch (c) {
            case '"':
                quoted += "\\\"";
                break;
            case '\\':
                quoted += "\\\\";
                break;
            case '\n':
                quoted += "\\n";
                break;
            case '\r':
                quoted += "\\r";
                break;
            case '\t':
                quoted += "\\t";
                break;
            default: {
                const auto byte = static_cast<unsigned char>(c);
                if (byte < 0x20) {
                    quoted += "\\u00";
                    quoted += kHex[byte >> 4U];
                    quoted += kHex[byte & 0xFU];
                } else {
                    quoted += c;
                }
            }
        }
    }
    return quoted + "\"";
}

std::string to_json(const core::Error& error) {
    return R"({"error":{"code":)" + json_string(core::to_string(error.code)) +
           ",\"message\":" + json_string(error.message) + "}}";
}

std::string groups_json() {
    return array(chem::kAllGroups,
                 [](chem::Group group) { return json_string(chem::to_string(group)); });
}

std::string url_json(std::string_view url) {
    return "{\"url\":" + json_string(url) + "}";
}

std::string to_json(std::span<const nmr1d::CarbonAssignment> peaks) {
    return array(peaks, [](const nmr1d::CarbonAssignment& peak) {
        return "{\"shift\":" + number(peak.shift) + ",\"types\":" + types(peak.possible) + "}";
    });
}

std::string to_json(const Suggestions& suggestions) {
    std::string text = "[";
    for (const auto group : chem::kAllGroups) {
        const auto items = group_evidence(suggestions, group);
        if (!items.empty()) {
            text += (text.size() > 1 ? "," : "") + std::string{R"({"group":)"} +
                    json_string(chem::to_string(group)) + R"(,"evidence":[)" + items + "]}";
        }
    }
    return text + "]";
}

std::string to_json(std::span<const io::ObservedPeak> peaks) {
    return array(peaks, [](const io::ObservedPeak& peak) {
        return "{\"shift\":" + number(peak.shift) +
               ",\"protons\":" + optional_count(filter::proton_count(peak)) + "}";
    });
}

std::string protons_json(const std::optional<std::vector<io::ObservedPeak>>& peaks) {
    return peaks.has_value() ? to_json(std::span<const io::ObservedPeak>{*peaks}) : "null";
}

std::string to_json(const Ranking& ranking) {
    const auto removal = [](const io::Removal& r) {
        return "{\"cid\":" + std::to_string(r.cid) + ",\"smiles\":" + json_string(r.smiles) +
               ",\"reason\":" + json_string(io::to_string(r.reason)) + "}";
    };
    return "{\"formula\":" + json_string(ranking.formula) + ",\"rows\":" + number(ranking.rows) +
           ",\"removed\":" + array(ranking.removed, removal) +
           ",\"candidates\":" + array(ranking.candidates, candidate) + "}";
}

}  // namespace arise::pipeline
