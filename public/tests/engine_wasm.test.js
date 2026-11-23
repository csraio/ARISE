// The real engine (build/web, from apps/wasm) called from JavaScript: values cross in both
// directions, answers parse, and every kind of bad input comes back as an EngineError.
// What the answers mean is tested in C++ (tests/pipeline), natively and in WebAssembly.
import assert from "node:assert/strict";
import { before, test } from "node:test";

import { EngineError } from "../js/engine.js";
import {
  benchmarkCase,
  benchmarkSearch,
  builtEngine,
  pastedSpectra,
  roundOneProtonText,
  roundOneProtonWindows,
} from "./fixtures.js";

let engine;
let benchmark;
before(async () => {
  engine = await builtEngine();
  benchmark = benchmarkCase();
});

const tierOrder = { active: 0, "back-burner": 1, eliminated: 2 };

function rankBenchmark(overrides = {}) {
  return engine.rank({
    formula: benchmark.formula,
    pubchemCsv: benchmarkSearch(benchmark),
    spectra: pastedSpectra(benchmark),
    present: [],
    absent: [],
    likely: [],
    useDept: true,
    ...overrides,
  });
}

function assertEngineError(call, code) {
  assert.throws(call, (error) => error instanceof EngineError && error.code === code);
}

test("searchUrl is a PubChem address naming the formula", () => {
  const url = engine.searchUrl(benchmark.formula);
  assert.ok(url.startsWith("https://pubchem.ncbi.nlm.nih.gov/"), url);
  assert.ok(url.includes(benchmark.formula), url);
});

test("groupNames lists distinct names", () => {
  const names = engine.groupNames();
  assert.ok(names.length > 0);
  assert.ok(names.every((name) => typeof name === "string" && name !== ""));
  assert.equal(new Set(names).size, names.length);
});

test("readSpectra classifies each benchmark peak as its type", () => {
  const peaks = engine.readSpectra(pastedSpectra(benchmark));
  assert.equal(peaks.length, benchmark.peaks.length);
  peaks.forEach((peak, i) => {
    assert.equal(peak.shift, Number(benchmark.peaks[i].shift));
    assert.deepEqual(peak.types, [benchmark.peaks[i].type]);
  });
});

test("readSpectra accepts CSV text and absent DEPT spectra", () => {
  const csv = "shift_ppm\n" + benchmark.peaks.map((p) => `${p.shift}\n`).join("");
  const peaks = engine.readSpectra({
    carbon: { format: "csv", text: csv },
    dept135: null,
    dept90: undefined,
  });
  assert.equal(peaks.length, benchmark.peaks.length);
  assert.ok(peaks.every((peak) => peak.types.length === 4));
});

test("suggest gives groups whose evidence lies inside its windows", () => {
  const groups = new Set(engine.groupNames());
  for (const useDept of [true, false]) {
    for (const suggestion of engine.suggest(pastedSpectra(benchmark), useDept)) {
      assert.ok(groups.has(suggestion.group), suggestion.group);
      assert.ok(suggestion.evidence.length > 0);
      for (const e of suggestion.evidence) {
        assert.ok(e.low <= e.shift && e.shift <= e.high, JSON.stringify(e));
      }
    }
  }
});

test("rank keeps every cleaned record, best tier first", () => {
  const ranking = rankBenchmark();
  assert.equal(ranking.formula, benchmark.formula);
  const dataRows = benchmarkSearch(benchmark).trimEnd().split("\n").length - 1;
  assert.equal(ranking.rows, dataRows);
  assert.equal(ranking.candidates.length + ranking.removed.length, ranking.rows);
  const tiers = ranking.candidates.map((c) => tierOrder[c.tier]);
  assert.ok(tiers.every((t) => t !== undefined));
  assert.ok(tiers.every((t, i) => i === 0 || tiers[i - 1] <= t), "tiers out of order");
  assert.ok(ranking.candidates.every((c) => Number.isInteger(c.cid) && c.cid > 0));
});

test("rank eliminates candidates without a group confirmed present", () => {
  const [first] = engine.suggest(pastedSpectra(benchmark), true);
  assert.ok(first, "the benchmark spectrum suggests no group");
  const ranking = rankBenchmark({ present: [first.group] });
  const eliminated = ranking.candidates.filter((c) => c.tier === "eliminated");
  assert.ok(eliminated.length > 0);
  assert.ok(eliminated.every((c) => c.reasons.length > 0));
});

test("the same call gives the same answer", () => {
  assert.deepEqual(rankBenchmark(), rankBenchmark());
});

test("engine errors keep their code", () => {
  assertEngineError(() => engine.searchUrl("not a formula"), "ParseError");
  const bad = pastedSpectra(benchmark);
  bad.carbon = { format: "pasted", text: "29,5\n" };
  assertEngineError(() => engine.readSpectra(bad), "ParseError");
  assertEngineError(() => rankBenchmark({ present: ["ketone"], absent: ["ketone"] }), "InvalidArgument");
});

test("values of the wrong kind are InvalidArgument, never a crash", () => {
  const spectra = pastedSpectra(benchmark);
  assertEngineError(() => engine.searchUrl(42), "InvalidArgument");
  assertEngineError(() => engine.searchUrl(undefined), "InvalidArgument");
  assertEngineError(() => engine.readSpectra(null), "InvalidArgument");
  assertEngineError(() => engine.readSpectra("text"), "InvalidArgument");
  assertEngineError(() => engine.readSpectra({ ...spectra, carbon: null }), "InvalidArgument");
  assertEngineError(
    () => engine.readSpectra({ ...spectra, carbon: { format: "xml", text: "" } }),
    "InvalidArgument",
  );
  assertEngineError(
    () => engine.readSpectra({ ...spectra, dept90: { format: "pasted", text: 5 } }),
    "InvalidArgument",
  );
  assertEngineError(() => engine.suggest(spectra, "yes"), "InvalidArgument");
  assertEngineError(() => rankBenchmark({ present: "ketone" }), "InvalidArgument");
  assertEngineError(() => rankBenchmark({ absent: ["no such group"] }), "InvalidArgument");
  assertEngineError(() => rankBenchmark({ absent: [7] }), "InvalidArgument");
  assertEngineError(() => rankBenchmark({ pubchemCsv: null }), "InvalidArgument");
  assertEngineError(() => rankBenchmark({ useDept: 1 }), "InvalidArgument");
  // The engine still answers after all of these.
  assert.ok(engine.groupNames().length > 0);
});

// ---- 1H ---------------------------------------------------------------------------------

function withProtons(text, format = "pasted") {
  return { ...pastedSpectra(benchmark), proton: { format, text } };
}

test("readProtons is null without a 1H spectrum and reads one when given", () => {
  assert.equal(engine.readProtons(pastedSpectra(benchmark)), null);
  const peaks = engine.readProtons(withProtons(roundOneProtonText()));
  const windows = roundOneProtonWindows();
  assert.equal(peaks.length, windows.length);
  peaks.forEach((peak, i) => assert.equal(peak.shift, Number(windows[i].median)));
  assert.ok(peaks.every((peak) => Number.isInteger(peak.protons)));
});

test("suggest gives 1H evidence for each round-1 window, with its precision", () => {
  const suggestions = engine.suggest(withProtons(roundOneProtonText()), true);
  for (const window of roundOneProtonWindows()) {
    const group = suggestions.find((s) => s.group === window.group);
    assert.ok(group, window.group);
    const evidence = group.evidence.filter((e) => e.spectrum === "1H");
    assert.ok(evidence.length > 0, window.group);
    assert.ok(evidence.every((e) => e.precision === window.precision && e.h_count === window.hCount));
  }
  // Without a 1H spectrum there is no 1H evidence.
  const carbonOnly = engine.suggest(pastedSpectra(benchmark), true);
  assert.ok(carbonOnly.every((s) => s.evidence.every((e) => e.spectrum === "13C")));
});

test("rank back-burners, never eliminates, candidates without a likely group", () => {
  const plain = rankBenchmark();
  const likely = rankBenchmark({ likely: ["methoxy"] });
  assert.equal(likely.candidates.length, plain.candidates.length);
  const eliminated = (r) => r.candidates.filter((c) => c.tier === "eliminated").length;
  assert.equal(eliminated(likely), eliminated(plain));
  const soft = likely.candidates.filter((c) =>
    c.reasons.includes("likely present (from 1H): methoxy; the candidate has none"),
  );
  assert.ok(soft.length > 0);
  assert.ok(soft.every((c) => c.tier === "back-burner"));
});

test("1H input errors are named", () => {
  assertEngineError(() => engine.readProtons(withProtons("3,8\n")), "ParseError");
  assert.throws(
    () => engine.readProtons(withProtons("3,8\n")),
    (error) => error.message.startsWith("1H spectrum"),
  );
  assertEngineError(
    () => engine.readProtons({ ...pastedSpectra(benchmark), proton: { format: "xml", text: "" } }),
    "InvalidArgument",
  );
  assertEngineError(() => rankBenchmark({ likely: "methoxy" }), "InvalidArgument");
  assertEngineError(() => rankBenchmark({ likely: ["methoxy"], absent: ["methoxy"] }), "InvalidArgument");
});

