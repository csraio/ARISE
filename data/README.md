# Data shipped with the engine

- `group_shifts.tsv`: each functional group's 13C hallmark shifts, measured over
  nmrshiftdb2's assigned spectra by `tools/measure_group_shifts` with ARISE's own group
  definitions. Compiled into the `filter` module at configure time; a test recomputes it
  from `tests/data/reference/nmrshiftdb2_13c_all.tsv`, so it cannot drift or be edited by
  hand.

Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made available here
under the nmrshiftdb2 Database License (`NMRSHIFTDB2_LICENSE.txt`). Clause 4.5 of that
licence requires software that relies on this data to be released under an OSI-approved
open-source licence; ARISE will be (user decision 2026-11-02).
