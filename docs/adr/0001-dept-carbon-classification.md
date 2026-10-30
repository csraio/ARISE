# ADR 0001: Classify 13C peaks by DEPT in a new `nmr1d` module

Date: 2026-10-30 · Status: accepted

## Context
Users want to know, for each 13C peak, whether the carbon is "primary, secondary,
tertiary or quaternary". DEPT measures something narrower: the number of hydrogens
attached to each carbon (C, CH, CH2, CH3). The common nicknames only match true carbon
degree for simple sp3 carbons. For example, an aromatic CH has two carbon neighbours,
and a ketone C=O, which DEPT reports as C, has two as well.

DEPT never shows carbons without hydrogens, so the broadband 13C spectrum is always
needed. DEPT-135 and DEPT-90 are each optional.

## Decision
- **Report attached-hydrogen counts** (`CarbonType`: C, CH, CH2, CH3), not carbon
  degree. Carbon degree needs a structure and belongs to `chem`, not to spectra.
- **Each peak gets a set of possible types** (`CarbonTypeSet`), narrowed by the spectra
  supplied, rather than a single forced answer:

  | Spectrum | Peak positive | Peak negative | Peak absent |
  |---|---|---|---|
  | DEPT-135 | {CH, CH3} | {CH2} | {C} |
  | DEPT-90 | {CH} | error | {C, CH2, CH3} |

  The supplied spectra are intersected. With 13C only, every type stays possible.
- **Matching:** each DEPT peak must match exactly one 13C peak within a tolerance
  (default 0.05 ppm, inclusive). Each 13C peak may match at most one peak per DEPT
  spectrum. Ambiguity, missing partners and contradictions are errors, never guesses.
- **New module `nmr1d`** (1D spectral interpretation) between `io` and `chem`. It needs
  the `io` peak types and nothing from RDKit. `io` stays parsing-only.
- **Peak lists:** 13C and DEPT use the shared peak-list CSV format. A DEPT list
  carries a signed `intensity` (the sign is the phase). A DEPT list may legitimately
  have no peaks; DEPT-90 of a molecule without CH is one example.

## Consequences
- `filter` can later compare the counts of C/CH/CH2/CH3 against each candidate
  structure, a strong and cheap constraint.
- Not detected: a DEPT-135 phased upside down when no DEPT-90 is supplied, and weak
  breakthrough peaks in DEPT-90 that users did not remove. A CDCl3 solvent peak is
  correctly reported as C, so users should leave solvent peaks out of the 13C list.
- The 0.05 ppm default is a convention, not a measured constant. Callers can lower it
  for crowded spectra.
