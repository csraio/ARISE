# Data shipped with the engine

- `solvent_peaks.tsv`: the residual-solvent and water 1H peaks of deuterated NMR solvents,
  read by `scripts/build_solvent_peaks.py` from the text layers of three published tables
  (Gottlieb et al., J. Org. Chem. 1997, 62, 7512; Fulmer et al., Organometallics 2010, 29,
  2176; the Cambridge Isotope Laboratories NMR Solvent Data Chart). A shift ships only where
  two tables agree within 0.01 ppm. Compiled into the `nmr1d` module at configure time; a
  test holds every shift to the values each table gives
  (`tests/data/reference/solvent_sources.tsv`). Only the numbers are kept, with citations,
  not the publications.
- `group_shifts.tsv`: each functional group's 13C hallmark shifts, measured over
  nmrshiftdb2's assigned spectra by `tools/measure_group_shifts` with ARISE's own group
  definitions. Compiled into the `filter` module at configure time; a test recomputes it
  from `tests/data/reference/nmrshiftdb2_13c_all.tsv`, so it cannot drift or be edited by
  hand.

`group_shifts.tsv` contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is
made available here under the nmrshiftdb2 Database License (`NMRSHIFTDB2_LICENSE.txt`).
Clause 4.5 of that licence requires software that relies on this data to be released under
an OSI-approved open-source licence; ARISE is released under Apache-2.0 (`LICENSE`).
