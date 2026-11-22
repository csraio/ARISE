#include "arise/pipeline/json.hpp"

#include <array>
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

std::string evidence(const filter::PeakEvidence& e) {
    return "{\"peak\":" + number(e.peak) + ",\"shift\":" + number(e.shift) +
           ",\"hallmark\":" + json_string(e.window.hallmark) +
           ",\"types\":" + types(e.window.types) + ",\"low\":" + number(e.window.low) +
           ",\"high\":" + number(e.window.high) + ",\"n\":" + number(e.window.n) +
           ",\"molecules\":" + number(e.window.molecules) +
           ",\"source\":" + json_string(source_name(e.window.source)) + "}";
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

std::string url_json(std::string_view url) {
    return "{\"url\":" + json_string(url) + "}";
}

std::string to_json(std::span<const nmr1d::CarbonAssignment> peaks) {
    return array(peaks, [](const nmr1d::CarbonAssignment& peak) {
        return "{\"shift\":" + number(peak.shift) + ",\"types\":" + types(peak.possible) + "}";
    });
}

std::string to_json(std::span<const filter::GroupSuggestion> suggestions) {
    return array(suggestions, [](const filter::GroupSuggestion& s) {
        return "{\"group\":" + json_string(chem::to_string(s.group)) +
               ",\"evidence\":" + array(s.evidence, evidence) + "}";
    });
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
