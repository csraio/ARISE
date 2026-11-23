// Shared by the Node tests: the repository's files and the built engine module.
import { readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

import { loadEngine } from "../js/engine.js";

export const repoRoot = resolve(dirname(fileURLToPath(import.meta.url)), "..", "..");

/** A repository file's text, by path from the repository root. */
export function repoText(path) {
  return readFileSync(resolve(repoRoot, path), "utf8");
}

/**
 * The engine built by `cmake --build --preset web`, or the module named by ARISE_ENGINE_JS.
 * Fails (never skips) when it is missing, so a test run cannot pass without the engine.
 */
export async function builtEngine() {
  const path =
    process.env.ARISE_ENGINE_JS ?? resolve(repoRoot, "build/web/apps/wasm/arise_engine.js");
  let factory;
  try {
    factory = (await import(pathToFileURL(path).href)).default;
  } catch (error) {
    throw new Error(
      `no engine at ${path}: build it with 'cmake --preset web && cmake --build --preset web' ` +
        `(after sourcing emsdk_env.sh) or set ARISE_ENGINE_JS (${error.message})`,
    );
  }
  return loadEngine(factory);
}

/**
 * The benchmark case (tests/data/reference/benchmark_spectrum.tsv, written by
 * scripts/build_reference_sets.py): its formula, and its peaks with one DEPT type each.
 */
export function benchmarkCase() {
  const peaks = [];
  let formula = "";
  for (const line of repoText("tests/data/reference/benchmark_spectrum.tsv").split("\n")) {
    if (line === "" || line.startsWith("#")) continue;
    const fields = line.split("\t");
    if (fields[0] === "formula") formula = fields[1];
    if (fields[0] === "peak") peaks.push({ shift: fields[1], type: fields[2] });
  }
  if (formula === "" || peaks.length === 0) throw new Error("benchmark case not found");
  return { formula, peaks };
}

/**
 * The spectra the case's types imply, pasted as "shift<TAB>intensity": DEPT-135 up for CH
 * and CH3, down for CH2, absent for C; DEPT-90 up for CH only.
 */
export function pastedSpectra(benchmark) {
  const carbon = benchmark.peaks.map((p) => `${p.shift}\n`).join("");
  const dept135 = benchmark.peaks
    .filter((p) => p.type !== "C")
    .map((p) => `${p.shift}\t${p.type === "CH2" ? -100 : 100}\n`)
    .join("");
  const dept90 = benchmark.peaks
    .filter((p) => p.type === "CH")
    .map((p) => `${p.shift}\t100\n`)
    .join("");
  return {
    carbon: { format: "pasted", text: carbon },
    dept135: { format: "pasted", text: dept135 },
    dept90: { format: "pasted", text: dept90 },
  };
}

/** The committed PubChem answer for the benchmark case's formula. */
export function benchmarkSearch(benchmark) {
  return repoText(`tests/data/pubchem_${benchmark.formula}.csv`);
}

/**
 * The round-1 1H windows of data/group_proton_shifts.tsv (written by
 * tools/measure_group_proton_shifts): [{group, hCount, median, precision}].
 */
export function roundOneProtonWindows() {
  const windows = [];
  for (const line of repoText("data/group_proton_shifts.tsv").split("\n")) {
    if (line === "" || line.startsWith("#")) continue;
    const f = line.split("\t");
    if (f[20] === "1") windows.push({ group: f[0], hCount: f[2], median: f[8], precision: Number(f[16]) });
  }
  if (windows.length === 0) throw new Error("no round-1 1H windows");
  return windows;
}

/** A pasted 1H list with one peak per round-1 window: its median shift, at its H count. */
export function roundOneProtonText() {
  const counts = { CH: 1, CH2: 2, CH3: 3 };
  return roundOneProtonWindows()
    .map((w) => `${w.median}\t${counts[w.hCount]}\n`)
    .join("");
}

