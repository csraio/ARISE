# ADR 0002: Expert-style pipeline: PubChem candidates, NMR pruning, engine + web client

Date: 2026-10-31 · Status: accepted

## Context
A chemist elucidates a structure by:
1. fixing the molecular formula;
2. collecting every known structure with that formula;
3. pruning that list with cheap 1D NMR evidence (counts of C/CH/CH2/CH3, carbonyl types,
   groups that must or must not be present);
4. only then spending effort on shift prediction for the dozen or so survivors.

ARISE should follow that order ("meta-computing"). Cheap, explainable evidence prunes first;
expensive physics runs only on the shortlist.

The original proposal had two candidate sources: a mass-indexed local library ("library
mode") and re-ranking of externally supplied lists. It also targeted an offline web page.

## Decision
1. **Candidates come from PubChem by molecular formula** (PUG REST `fastformula`,
   `ConnectivitySMILES`). PubChem is free, needs no account, and has about 120M
   compounds. MolPort is commercial and lists only purchasable compounds. Imported
   `.smi` lists remain supported. Built-in constrained isomer generation is deferred,
   for compounds PubChem lacks.
2. **Two deliverables with a clean boundary.**
   - **Engine:** portable C++20 (native and wasm). Parsing, cleaning, pruning and later
     the chemistry. Text in, structured results out; no network, no files, no browser
     code.
   - **Web client:** static files, with modern JavaScript (ES modules) for the UI, the
     PubChem `fetch()` and file upload/download. No backend. It is served from the
     user's machine or hosted publicly, with identical files either way. The native CLI
     uses saved files and has no network code.
3. **No new dependencies for the pruning stage.** The engine reads PubChem's Kekulé
   SMILES with a small in-house reader (`chem`), instead of wrapping RDKit. Whether RDKit
   is needed for 3D conformers is decided at the anisotropy stage.
4. **Pruning evidence:** C/CH/CH2/CH3 counts (from `nmr1d::classify_carbons`),
   symmetry-unique carbons versus 13C peaks, and functional groups. ARISE suggests
   groups from shifts with their evidence; only groups the chemist confirms become hard
   filters. Every elimination states its reason.
5. **Baseline shift prediction stays an owned HOSE-code predictor** (not nmrdb.org,
   whose engine is the proprietary SPINUS network), so the 3D anisotropy correction is
   defined against a known baseline. Its design is on hold until the user provides
   context.

## Evidence (live requests, 2026-10-31)
- PubChem answers browser requests from any origin (`access-control-allow-origin: *`),
  so no user needs to relax browser security.
- One request per formula returns every match:

  | Formula | Time | Rows | Distinct |
  |---|---|---|---|
  | C4H8O2 | 1.3 s | 494 | 211 |
  | C8H8O2 | 1.5 s | 1,239 | 925 |
  | C10H12O2 | 2.8 s | 5,544 | 3,986 |
  | C15H24O | 4.1 s | 8,833 | 5,706 |

- Raw results contain duplicates (butanoic acid appears 38 times for C4H8O2),
  multi-component entries, charged species and radicals. Cleaning is the first
  pruning stage.

## Consequences
- **Speed target:** formula to pruned list in about 2–5 s, dominated by PubChem. The
  engine must clean and prune 10,000 candidates in under 100 ms (benchmark). Huge
  result sets fall back to fetching CIDs first, then SMILES in POST batches of about
  1,000, at most 5 requests per second.
- **Determinism is defined for a given candidate list.** PubChem changes over time, so
  results record the formula, query date and retrieved list, which can be re-imported.
- **Compounds absent from PubChem** (new natural products, new syntheses) cannot be
  found until built-in generation exists.
- **Superseded parts of the proposal:** library mode, the offline web page, RDKit in the
  pruning stages, and the "re-rank only" scope. The proposal document needs a revision.
