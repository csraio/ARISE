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

- `hybridization_shifts.tsv`: the 13C shifts of sp, sp2 and sp3 carbons, measured over
  nmrshiftdb2's assigned spectra by `tools/measure_hybridization_shifts` with ARISE's own
  hybridization rule (`chem::carbon_hybridization`), for suggesting sp2/sp3 peak counts
  when there is no DEPT. A test recomputes it from
  `tests/data/reference/nmrshiftdb2_13c_all.tsv`, and another holds it to the builder's
  independent RDKit measurement (`tests/data/reference/measurements.tsv`).

- `group_proton_shifts.tsv`: each functional group's 1H windows (the H on its key carbons,
  split by the carbon's number of H), with each window's measured precision and the group's
  base rate, measured over nmrshiftdb2's assigned 1H spectra by
  `tools/measure_group_proton_shifts` with ARISE's own group definitions. Its `round` column
  marks the windows the chemist enabled (2026-11-22: methoxy and benzene ring). A test
  recomputes it from `tests/data/reference/nmrshiftdb2_1h_all.tsv`, and another holds the
  round-1 windows to the builder's independent RDKit measurement.

`group_shifts.tsv`, `group_proton_shifts.tsv` and `hybridization_shifts.tsv` contain
information from nmrshiftdb2 (www.nmrshiftdb.org), which is
made available here under the nmrshiftdb2 Database License (`NMRSHIFTDB2_LICENSE.txt`).
Clause 4.5 of that licence requires software that relies on this data to be released under
an OSI-approved open-source licence; ARISE is released under Apache-2.0 (`LICENSE`).
