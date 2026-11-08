# ADR 0003: Separate `hsqc` and `hmbc` modules, before `chem`

Date: 2026-11-05 · Status: accepted

## Context
The module chain in CLAUDE.md had one `nmr2d` module near the end, after `aniso`, for
"HSQC/HMBC scoring". Two things changed that:

- HSQC is needed early, in 1D-style pruning, not only for final scoring. A 1H peak
  with no HSQC cross peak, which is also not a solvent or water peak, is suggested as
  exchangeable (H on O, N, S or Se). A count of such protons that the chemist confirms is
  a hard fact about the candidates (plan 3x).
- Reading a spectrum and comparing it with a candidate are different jobs. Deciding
  which 1H peaks have an HSQC partner needs only peak lists. Checking a candidate's C–H
  pairs needs a molecule.

The user decided on 2026-11-01 that HSQC and HMBC become their own modules.

## Decision
- **Two sibling modules after `nmr1d` and before `chem`:**
  - `hsqc`: HSQC peak sets and 1H↔HSQC matching (`classify_protons`,
    `suggest_exchangeable_count`).
  - `hmbc`: planned only. It gets code, and a library, when the HMBC stage starts.
  - Neither depends on the other. Both work on spectra alone, without molecules.
- **`nmr2d` is removed** from the chain.
- **Candidate-vs-spectrum scoring needs molecules,** so it lives in `filter` (or a later
  scoring module). That code calls `hsqc` and `hmbc`.

New chain: `core` → `io` → `nmr1d` → {`hsqc`, `hmbc`} → `chem` → `filter` → `predict`
→ `aniso` → `pipeline` → `apps/*`.

## Consequences
- `hsqc` can use `nmr1d`'s solvent table. `chem` and everything after it can use
  `hsqc`.
- The first `hsqc` code is the exchangeable-proton suggestion. The 1H↔HSQC matching
  tolerance is a convention the chemist set: 0.15 ppm by default, adjustable. It was
  measured on BMRB experimental spectra (`tests/data/reference/bmrb_*.tsv`).
- No empty `hmbc` library is added. CLAUDE.md lists the module so its place is fixed.
