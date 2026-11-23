// The engine as the web client calls it (Emscripten embind). Every function takes plain
// JavaScript values and returns JSON text (arise/pipeline/json.hpp); a failure comes back
// as {"error":{"code":...,"message":...}}, never as an exception or an abort, including
// for a value of the wrong JavaScript type. The work itself is arise::pipeline; this file
// only converts values.
//
// A spectrum is {format: "csv" | "pasted", text: string}; the spectra are
// {carbon: spectrum, dept135: spectrum | null, dept90: spectrum | null,
//  proton: spectrum | null} (proton may also be left out).

#include "arise/core/error.hpp"
#include "arise/core/result.hpp"
#include "arise/pipeline/json.hpp"
#include "arise/pipeline/pipeline.hpp"

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <optional>
#include <string>
#include <utility>

namespace {

using arise::core::Error;
using arise::core::ErrorCode;
using arise::core::Result;
using arise::pipeline::SpectraText;
using arise::pipeline::SpectrumText;
using emscripten::val;

Error invalid(const std::string& message) {
    return {.code = ErrorCode::InvalidArgument, .message = message};
}

// Whether `value` is a JavaScript object (not null).
bool is_object(const val& value) {
    return !value.isNull() && value.typeOf().as<std::string>() == "object";
}

Result<std::string> string_of(const val& value, const std::string& what) {
    if (!value.isString()) {
        return invalid(what + " must be text");
    }
    return value.as<std::string>();
}

Result<SpectrumText> spectrum_of(const val& value, const std::string& what) {
    if (!is_object(value)) {
        return invalid(what + " must be {format, text}");
    }
    const auto format_name = string_of(value["format"], what + ".format");
    if (!format_name.ok()) {
        return format_name.error();
    }
    const auto format = arise::pipeline::text_format_named(format_name.value());
    if (!format.has_value()) {
        return invalid(what + ".format must be \"csv\" or \"pasted\", not \"" +
                       format_name.value() + "\"");
    }
    auto text = string_of(value["text"], what + ".text");
    if (!text.ok()) {
        return std::move(text).error();
    }
    return SpectrumText{.format = *format, .text = std::move(text).value()};
}

Result<std::optional<SpectrumText>> optional_spectrum_of(const val& value,
                                                         const std::string& what) {
    if (value.isNull() || value.isUndefined()) {
        return std::optional<SpectrumText>{};
    }
    auto spectrum = spectrum_of(value, what);
    if (!spectrum.ok()) {
        return std::move(spectrum).error();
    }
    return std::optional{std::move(spectrum).value()};
}

Result<SpectraText> spectra_of(const val& value) {
    if (!is_object(value)) {
        return invalid("spectra must be {carbon, dept135, dept90}");
    }
    auto carbon = spectrum_of(value["carbon"], "spectra.carbon");
    if (!carbon.ok()) {
        return std::move(carbon).error();
    }
    auto dept135 = optional_spectrum_of(value["dept135"], "spectra.dept135");
    if (!dept135.ok()) {
        return std::move(dept135).error();
    }
    auto dept90 = optional_spectrum_of(value["dept90"], "spectra.dept90");
    if (!dept90.ok()) {
        return std::move(dept90).error();
    }
    auto proton = optional_spectrum_of(value["proton"], "spectra.proton");
    if (!proton.ok()) {
        return std::move(proton).error();
    }
    return SpectraText{.carbon = std::move(carbon).value(),
                       .dept135 = std::move(dept135).value(),
                       .dept90 = std::move(dept90).value(),
                       .proton = std::move(proton).value()};
}

// The groups named in the JavaScript array `value`.
Result<std::vector<arise::chem::Group>> groups_of(const val& value, const std::string& what) {
    if (!val::global("Array").call<bool>("isArray", value)) {
        return invalid(what + " must be a list of group names");
    }
    std::vector<arise::chem::Group> groups;
    const auto length = value["length"].as<unsigned>();
    for (unsigned i = 0; i < length; ++i) {
        const auto name = string_of(value[i], what + " entry");
        if (!name.ok()) {
            return name.error();
        }
        const auto group = arise::pipeline::group_named(name.value());
        if (!group.has_value()) {
            return invalid(what + ": no group named \"" + name.value() + "\"");
        }
        groups.push_back(*group);
    }
    return groups;
}

Result<bool> flag_of(const val& value, const std::string& what) {
    if (!value.isTrue() && !value.isFalse()) {
        return invalid(what + " must be true or false");
    }
    return value.isTrue();
}

// `result` as JSON, written by `write`, or its error.
template <typename T, typename Write>
std::string json_of(const Result<T>& result, Write write) {
    return result.ok() ? write(result.value()) : arise::pipeline::to_json(result.error());
}

std::string search_url(const val& formula) {
    const auto text = string_of(formula, "formula");
    if (!text.ok()) {
        return arise::pipeline::to_json(text.error());
    }
    return json_of(arise::pipeline::search_url(text.value()),
                   [](const std::string& url) { return arise::pipeline::url_json(url); });
}

std::string group_names() {
    return arise::pipeline::groups_json();
}

std::string read_spectra(const val& spectra) {
    const auto parsed = spectra_of(spectra);
    if (!parsed.ok()) {
        return arise::pipeline::to_json(parsed.error());
    }
    return json_of(arise::pipeline::read_spectra(parsed.value()), [](const auto& peaks) {
        return arise::pipeline::to_json(std::span<const arise::nmr1d::CarbonAssignment>{peaks});
    });
}

std::string suggest(const val& spectra, const val& use_dept) {
    const auto parsed = spectra_of(spectra);
    const auto dept = flag_of(use_dept, "useDept");
    if (!parsed.ok() || !dept.ok()) {
        return arise::pipeline::to_json(parsed.ok() ? dept.error() : parsed.error());
    }
    return json_of(arise::pipeline::suggest(parsed.value(), {.use_dept = dept.value()}, {}),
                   [](const arise::pipeline::Suggestions& suggestions) {
                       return arise::pipeline::to_json(suggestions);
                   });
}

std::string read_protons(const val& spectra) {
    const auto parsed = spectra_of(spectra);
    if (!parsed.ok()) {
        return arise::pipeline::to_json(parsed.error());
    }
    return json_of(arise::pipeline::read_protons(parsed.value()),
                   [](const auto& peaks) { return arise::pipeline::protons_json(peaks); });
}

// The confirmed groups from the JavaScript arrays `present`, `absent` and `likely` (groups
// confirmed from 1H evidence only: soft, filter::ConfirmedGroups::likely_present).
Result<arise::filter::ConfirmedGroups> confirmed_of(const val& present, const val& absent,
                                                    const val& likely) {
    auto yes = groups_of(present, "present");
    if (!yes.ok()) {
        return std::move(yes).error();
    }
    auto no = groups_of(absent, "absent");
    if (!no.ok()) {
        return std::move(no).error();
    }
    auto maybe = groups_of(likely, "likely");
    if (!maybe.ok()) {
        return std::move(maybe).error();
    }
    return arise::filter::ConfirmedGroups{.present = std::move(yes).value(),
                                          .absent = std::move(no).value(),
                                          .likely_present = std::move(maybe).value()};
}

std::string rank(const val& formula, const val& pubchem_csv, const val& spectra, const val& present,
                 const val& absent, const val& likely, const val& use_dept) {
    const auto formula_text = string_of(formula, "formula");
    const auto csv = string_of(pubchem_csv, "pubchemCsv");
    const auto parsed = spectra_of(spectra);
    const auto confirmed = confirmed_of(present, absent, likely);
    const auto dept = flag_of(use_dept, "useDept");
    for (const Error* error :
         {formula_text.ok() ? nullptr : &formula_text.error(), csv.ok() ? nullptr : &csv.error(),
          parsed.ok() ? nullptr : &parsed.error(), confirmed.ok() ? nullptr : &confirmed.error(),
          dept.ok() ? nullptr : &dept.error()}) {
        if (error != nullptr) {
            return arise::pipeline::to_json(*error);
        }
    }
    return json_of(
        arise::pipeline::rank(formula_text.value(), csv.value(), parsed.value(), confirmed.value(),
                              {.use_dept = dept.value()}),
        [](const arise::pipeline::Ranking& ranking) { return arise::pipeline::to_json(ranking); });
}

}  // namespace

EMSCRIPTEN_BINDINGS(arise_engine) {
    emscripten::function("searchUrl", &search_url);
    emscripten::function("groupNames", &group_names);
    emscripten::function("readSpectra", &read_spectra);
    emscripten::function("suggest", &suggest);
    emscripten::function("readProtons", &read_protons);
    emscripten::function("rank", &rank);
}
