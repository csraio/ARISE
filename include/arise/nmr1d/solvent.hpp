#pragma once

#include "arise/core/result.hpp"
#include "arise/core/units.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace arise::nmr1d {

/// The deuterated NMR solvents whose residual 1H peaks at least two published tables agree
/// on (data/solvent_peaks.tsv, read from Gottlieb 1997, Fulmer 2010 and the CIL chart by
/// scripts/build_solvent_peaks.py).
enum class Solvent : std::uint8_t {
    ThfD8,
    CD2Cl2,
    CDCl3,
    TolueneD8,
    C6D6,
    AcetoneD6,
    DmsoD6,
    CD3CN,
    TfeD3,
    CD3OD,
    D2O,
};

/// Every Solvent, in the column order of Fulmer et al. 2010, Table 1 (for menus).
inline constexpr std::array kAllSolvents{
    Solvent::ThfD8, Solvent::CD2Cl2,    Solvent::CDCl3,  Solvent::TolueneD8,
    Solvent::C6D6,  Solvent::AcetoneD6, Solvent::DmsoD6, Solvent::CD3CN,
    Solvent::TfeD3, Solvent::CD3OD,     Solvent::D2O,
};

/// The solvent's name as the published tables print it ("CDCl3", "(CD3)2SO", "THF-d8").
[[nodiscard]] std::string_view to_string(Solvent solvent) noexcept;

/// Where a solvent's own signals fall in a 1H spectrum.
struct SolventPeaks {
    Solvent solvent = Solvent::CDCl3;
    /// The PubChem compound the tables' heading was identified with.
    std::size_t pubchem_cid = 0;
    /// Residual (partly protonated solvent) peaks, ascending; never empty.
    std::vector<core::Ppm> residual;
    /// Water (H2O, and HDO where the tables give it), ascending; empty when the tables do
    /// not agree on a value. Water moves with temperature and concentration, so a peak here
    /// is "possibly water", never certainly.
    std::vector<core::Ppm> water;
    /// True when the solvent carries deuterium on O, N or S (D2O, CD3OD, TFE-d3): it swaps a
    /// solute's OH/NH/SH protons for deuterium, so they may not be seen at all.
    bool exchanges_labile_h = false;

    friend bool operator==(const SolventPeaks&, const SolventPeaks&) = default;
};

/// Reads a table in the format of data/solvent_peaks.tsv: '#' comment lines, then
/// tab-separated rows of name (as to_string), PubChem CID, residual shifts and water shifts
/// (each comma-separated, ascending; water may be '-' for none), and "yes" or "no" for
/// exchanges_labile_h. Rows are returned in file order.
///
/// Errors (ParseError, citing the line): a row without 5 fields, an unknown solvent name, a
/// solvent listed twice, a CID that is not a positive whole number, no residual shift, a
/// shift that is not a plain decimal, outside [io::kMinProtonShift, io::kMaxProtonShift] or
/// out of ascending order, an exchange field other than "yes"/"no", and a table with no rows.
[[nodiscard]] core::Result<std::vector<SolventPeaks>> parse_solvent_peaks(std::string_view text);

/// The peaks ARISE ships for one solvent (data/solvent_peaks.tsv, compiled in at configure
/// time and read with parse_solvent_peaks). NotFound if the table lacks the solvent.
[[nodiscard]] core::Result<SolventPeaks> solvent_peaks(Solvent solvent);

}  // namespace arise::nmr1d
