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
