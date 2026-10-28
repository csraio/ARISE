# ARISE — Anisotropy-Resolved Isomer Scoring Engine

Deterministic, explainable NMR structure-ranking engine. It re-ranks candidate
small molecules (mostly regioisomers) using 1H shifts, a HOSE-code baseline, and
a 3D ring-current/anisotropy correction, plus optional HSQC/HMBC scoring.
Ships as a native CLI and as a static WebAssembly page.
Design and scope: `docs/ARISE_Project_Proposal.pdf` (source: `.html`).

## Golden rules
1. The engine (`include/`, `src/`) is idiomatic, portable C++20. It must compile
   natively and with Emscripten without code changes.
2. Every module has GoogleTest tests. No commit without a green build and
   passing tests (enforced by `.githooks/pre-commit`; never use `--no-verify`).
3. Python is allowed only in `research/` (throwaway experiments) and `scripts/`
   (development-time data builders). The engine never calls, embeds or needs Python;
   the CMake build must work on a machine without Python.
4. Same input gives the same output, byte-for-byte, on native and WebAssembly.

## Repository layout
- `include/arise/<module>/` public headers · `src/<module>/` implementation
- `tests/<module>/` GoogleTest files mirroring `src/` · `tests/golden/` reference data
- `bench/` Google Benchmark · `apps/cli/` native CLI · `apps/wasm/` Emscripten bindings
- `public/` static web page · `scripts/` Python data builders · `research/` Python experiments
- `third-party/` pinned dependencies (do not edit) · `docs/adr/` design decisions

Modules (namespace `arise::<module>`), dependencies point downward only:
`core` (types, units, Result) → `io` (CSV, packed binary formats) → `chem` (RDKit wrapper)
→ `filter` (mass, SMARTS) → `predict` (HOSE 1H) → `aniso` (conformers, ring currents)
→ `nmr2d` (HSQC/HMBC scoring) → `pipeline` → `apps/*`.

## Build & test
    git config core.hooksPath .githooks      # once per clone
    cmake --preset dev && cmake --build --preset dev
    ctest --preset dev --output-on-failure
    cmake --preset asan ...                  # sanitizers
    emcmake cmake --preset wasm ...          # WebAssembly build

## C++ conventions
- C++20, no compiler extensions. Allowed: concepts, `std::span`, `<bit>`,
  designated initializers, `constexpr`. Avoid for now: modules, coroutines,
  `std::format` (Apple clang 14 lacks it; use fmt), and `<ranges>` beyond the basics.
- RAII only; no raw `new`/`delete`; ownership with `std::unique_ptr`/values.
  Pass views (`std::span`, `std::string_view`) in, return values out.
- `const` by default, `[[nodiscard]]` on functions returning results, `enum class`,
  `noexcept` where true. No macros except include guards (`#pragma once`).
- Errors: `arise::Result<T>` (value or `Error`). RDKit exceptions are caught in
  `chem` and converted; exceptions never cross ARISE module boundaries.
- Strong unit types: `Ppm`, `Hz`, `Dalton`, `Angstrom`. Never pass bare doubles
  for physical quantities across module boundaries.
- Naming: `snake_case` functions/variables, `PascalCase` types,
  `kConstant` constants, `member_` private members, files `snake_case.hpp/.cpp`.
- Formatting by `.clang-format`; lint by `.clang-tidy`. Do not hand-format.

## Portability (WebAssembly) rules
- Core modules: no file I/O, no threads, no globals with mutable state, no
  platform headers, no x86 intrinsics/inline asm, no `long double`.
- Don't assume `sizeof(size_t) == 8` (wasm32). Binary formats use fixed-width
  little-endian fields with a version header.
- Determinism: seeded RNGs (seed from config), stable sorts with explicit
  tie-breaks, never let `unordered_*` iteration order reach output.
- No `-ffast-math`. Floating-point tests compare with explicit tolerances.

## Testing
- Each new public function gets tests in `tests/<module>/` in the same commit.
- Name tests `TEST(ModuleName, BehaviorUnderCondition)`.
- Physics and chemistry code: test against hand-computed values and known
  literature cases; cite the source in a comment.
- Golden tests compare the pipeline to `tests/golden/` within stated tolerances.
  Regenerate golden files only on purpose, with a commit message explaining why.
- Speed targets (Stages 1–2: 500 candidates < 50 ms; full: 12 candidates < 1 s)
  are guarded by benchmarks in `bench/`.

## Data files
- Built by `scripts/*.py`, consumed by C++. The C++ side owns and validates
  the format; each file records its generator script and source-data version.
- Large generated data is not committed; document how to rebuild it.

## Dependencies
- Pinned exact versions via CMake FetchContent in `third-party/`: RDKit (only
  GraphMol, SmilesParse, Substruct, DistGeomHelpers, ForceFieldHelpers; Python
  wrappers OFF), GoogleTest, Google Benchmark, fmt. Ask before adding any new one.

## Git workflow
- Branch from `main`; small commits; Conventional Commits
  (`feat(aniso): add Johnson–Bovey ring current`).
- Never commit build output, data dumps or secrets. Never bypass hooks.
- Significant design choices get a short ADR in `docs/adr/`.

## For Claude
- Tests come first. Before implementing or changing code, write rigorous
  GoogleTest cases for every module touched: normal behavior, edge cases,
  invalid input, and known reference values. Then implement against them.
- In auto mode, never ship or run a full build until those tests exist and
  cover the change. A change "works" only when its tests prove it, not
  because a build succeeded.
- Report test failures with their output; never weaken or delete a test to
  make it pass.
- Keep changes inside one module where possible; follow existing patterns.
- When the science is uncertain (shift models, constants), say so and cite sources
  rather than guessing values.
