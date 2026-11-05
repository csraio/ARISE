# ARISE — Anisotropy-Resolved Isomer Scoring Engine

Deterministic, explainable NMR structure-elucidation engine that works like an expert:
candidates with the known molecular formula come from PubChem, cheap 1D NMR evidence
(DEPT counts, symmetry, confirmed functional groups) prunes them, and only the
survivors get shift prediction (HOSE baseline + 3D ring-current/anisotropy
correction) and optional HSQC/HMBC scoring.
Two deliverables: the portable C++ engine (native CLI and WebAssembly) and a static
web client (JavaScript UI). Design: `docs/adr/0002-expert-pipeline-pubchem-web.md`;
original proposal `docs/ARISE_Project_Proposal.pdf` (partly superseded by ADR 0002).

## Golden rules
1. The engine (`include/`, `src/`) is idiomatic, portable C++20. It must compile
   natively and with Emscripten without code changes.
2. Every module has GoogleTest tests. No commit without a green build and
   passing tests (enforced by `.githooks/pre-commit`; never use `--no-verify`).
3. Python is allowed only in `research/` (throwaway experiments) and `scripts/`
   (development-time data builders). The engine never calls, embeds or needs Python;
   the CMake build must work on a machine without Python.
4. Same input gives the same output, byte-for-byte, on native and WebAssembly. Input
   includes the candidate list: PubChem changes over time, so results record the
   formula, query date and retrieved list.
5. The engine never touches the network or files. Only the app layer does: the web
   client fetches PubChem and handles uploads/downloads, the CLI reads saved files.
   The engine takes text in and returns results.

## Repository layout
- `include/arise/<module>/` public headers · `src/<module>/` implementation
- `tests/<module>/` GoogleTest files mirroring `src/` · `tests/golden/` reference data
- `bench/` Google Benchmark · `apps/cli/` native CLI · `apps/wasm/` Emscripten bindings
- `public/` web client (static HTML + JavaScript ES modules, no backend)
- `scripts/` Python data builders · `research/` Python experiments
- `data/` data compiled into the engine at configure time (e.g. measured group shifts),
  each file generated and test-checked, with its source licence beside it
- `tools/` development-only C++ tools (e.g. measuring group shifts with ARISE's own
  definitions); built with the tests, never part of the engine, CLI or web app
- `third-party/` pinned dependencies (do not edit) · `docs/adr/` design decisions

Modules (namespace `arise::<module>`), dependencies point downward only:
`core` (types, units, Result) → `io` (CSV, packed binary formats)
→ `nmr1d` (1D spectral interpretation, e.g. DEPT, solvents)
→ {`hsqc` (1H↔HSQC matching, exchangeable protons), `hmbc` (planned)} — siblings,
neither depends on the other → `chem` (in-house SMILES reader, molecule graph)
→ `filter` (DEPT counts, symmetry, functional groups; candidate-vs-spectrum checks)
→ `predict` (HOSE baseline) → `aniso` (conformers, ring currents) → `pipeline`
→ `apps/*`. See `docs/adr/0003-hsqc-hmbc-modules.md`.

## Build & test
    git config core.hooksPath .githooks      # once per clone
    sh scripts/install_dev_tools.sh          # once per machine: pinned clang-format/clang-tidy
    cmake --preset dev && cmake --build --preset dev
    ctest --preset dev --output-on-failure
    cmake --preset asan ...                  # sanitizers
    sh scripts/install_emsdk.sh              # once per machine: pinned Emscripten + Node 22
    . ~/.local/share/emsdk/emsdk_env.sh      # per shell, then the wasm preset:
    cmake --preset wasm && cmake --build --preset wasm && ctest --preset wasm
The wasm preset runs the tests under Node with UBSan, SAFE_HEAP, ASSERTIONS=2 and
stack-overflow checks. CI runs dev, asan and release on Linux and macOS, plus wasm.

## C++ conventions
- C++20, no compiler extensions. Allowed: concepts, `std::span`, `<bit>`,
  designated initializers, `constexpr`. Avoid for now: modules, coroutines,
  `std::format` (Apple clang 14 lacks it; use fmt), and `<ranges>` beyond the basics.
- RAII only; no raw `new`/`delete`; ownership with `std::unique_ptr`/values.
  Pass views (`std::span`, `std::string_view`) in, return values out.
- Index with `operator[]`, not `.at()` (which throws). Debug builds harden the standard
  library (`_GLIBCXX_ASSERTIONS`, `_LIBCPP_HARDENING_MODE`), so out-of-range access traps
  in every test run instead of returning garbage.
- `const` by default, `[[nodiscard]]` on functions returning results, `enum class`,
  `noexcept` where true. No macros except include guards (`#pragma once`).
- Errors: `arise::Result<T>` (value or `Error`). Exceptions never cross ARISE module
  boundaries; any third-party exception is caught and converted where it is called.
- Strong unit types: `Ppm`, `Hz`, `Dalton`, `Angstrom`. Never pass bare doubles
  for physical quantities across module boundaries.
- Naming: `snake_case` functions/variables, `PascalCase` types,
  `kConstant` constants, `member_` private members, files `snake_case.hpp/.cpp`.
- Formatting by `.clang-format`; lint by `.clang-tidy` (`scripts/run_clang_tidy.sh`).
  Do not hand-format. Tool versions are pinned in `scripts/dev-tools-requirements.txt`;
  CI checks both.

## Portability (WebAssembly) rules
- Core modules: no file I/O, no threads, no globals with mutable state, no
  platform headers, no x86 intrinsics/inline asm, no `long double`.
- Don't assume `sizeof(size_t) == 8` (wasm32). Binary formats use fixed-width
  little-endian fields with a version header.
- Determinism: seeded RNGs (seed from config), stable sorts with explicit
  tie-breaks, never let `unordered_*` iteration order reach output.
- No `-ffast-math`, and no FMA contraction (`-ffp-contract=off` is set globally;
  arm64 clang would otherwise fuse `a*b+c` and change last bits). Floating-point
  tests compare with explicit tolerances; `tests/core/platform_test.cpp` guards
  the IEEE-754 assumptions.

## Testing
- Each new public function gets tests in `tests/<module>/` in the same commit.
- Name tests `TEST(ModuleName, BehaviorUnderCondition)`.
- Physics and chemistry code: test against hand-computed values and known
  literature cases; cite the source in a comment.
- Golden tests compare the pipeline to `tests/golden/` within stated tolerances.
  Regenerate golden files only on purpose, with a commit message explaining why.
- Speed targets, guarded by benchmarks in `bench/`: cleaning + pruning 10,000
  candidates < 100 ms (native and wasm); full chemistry on 12 candidates < 1 s.
  End to end, formula to pruned list should take about 2–5 s, mostly PubChem.
- The commit gate tests exactly the staged snapshot (exported to `.git/arise-gate/`),
  so unstaged edits and untracked files can neither break nor rescue a commit.
  The gate itself is tested by `scripts/verify_commit_gate.sh` (also in CI).
- Mutation testing (Mull) is planned for CI once more modules exist.

## Data files
- Built by `scripts/*.py`, consumed by C++. The C++ side owns and validates
  the format; each file records its generator script and source-data version.
- Large generated data is not committed; document how to rebuild it.

## Dependencies
- As few as possible. Pinned exact versions via CMake FetchContent in `third-party/`:
  GoogleTest (tests only); Google Benchmark planned (bench only). Pinned tools:
  clang-format/clang-tidy, Emscripten. No RDKit, libcurl or JS framework; whether
  RDKit is needed for 3D conformers is decided at the anisotropy stage.
  Ask before adding any new dependency.

## Git workflow
- Branch from `main`; small commits; Conventional Commits
  (`feat(aniso): add Johnson–Bovey ring current`).
- Never commit build output, data dumps or secrets. Never bypass hooks.
- Significant design choices get a short ADR in `docs/adr/`.
- CI (`.github/workflows/ci.yml`) builds and tests every push and PR on Linux and
  macOS (`dev`, `asan`, `release`), plus wasm, lint and the commit-gate self-test.
  `main` only takes changes through green PRs.

## Safe experiments & shell hygiene
Lessons from a real incident: a deliberately broken test was "reverted" with
`git checkout` on a file that had never been committed, so nothing was restored,
and the failure was hidden by `;` chaining and `2>/dev/null`.

- **Isolate experiments.** Never edit, break or stage live working-tree files to
  test tooling (hooks, CI, test harnesses). Use `scripts/verify_commit_gate.sh` or
  a throwaway `git worktree add --detach "$(mktemp -d)" HEAD`.
- **Commit or stash before experimenting**, so every file has a git copy to
  return to. `git checkout -- <file>` cannot restore a file git has never seen.
- **Chain dependent steps with `&&`**, not `;`, and use `set -eu` in scripts,
  so the first failure stops everything after it.
- **Never silence cleanup or revert commands.** No `2>/dev/null`, `|| true` or
  `>/dev/null 2>&1` on anything whose failure would leave the tree changed.
- **Always clean up.** Scripts that create temp dirs, worktrees or branches remove
  them in a `trap '…' EXIT`.
- **Verify state, not strings.** After a revert or experiment, prove the tree is
  unchanged with `git status --porcelain` (empty) and `git diff --exit-code`,
  not by grepping for a pattern.
- **Scripts are POSIX `sh`; test them with `dash`.** On macOS `sh` is really bash, so
  bashisms pass locally and fail in CI (Ubuntu's `/bin/sh` is dash). Run
  `dash scripts/<name>.sh` before committing a script. `sh` has no `pipefail`:
  capture a command's output before piping it, so a crash cannot hide behind `| head`.
- **Stage files explicitly.** Prefer `git add <paths>` over `git add -A`, and
  check `git status --short` before every commit.

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
