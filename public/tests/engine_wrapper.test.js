// public/js/engine.js with a stand-in engine module: answers are parsed, and error answers
// become EngineError. The real engine is tested in engine_wasm.test.js.
import assert from "node:assert/strict";
import { test } from "node:test";

import { EngineError, parseAnswer, wrapEngine } from "../js/engine.js";

test("parseAnswer returns the parsed value", () => {
  assert.deepEqual(parseAnswer('{"url":"u"}'), { url: "u" });
  assert.deepEqual(parseAnswer('["a","b"]'), ["a", "b"]);
  assert.deepEqual(parseAnswer("[]"), []);
});

test("parseAnswer throws EngineError for an error answer", () => {
  assert.throws(
    () => parseAnswer('{"error":{"code":"ParseError","message":"line 2: bad"}}'),
    (error) =>
      error instanceof EngineError &&
      error.code === "ParseError" &&
      error.message === "line 2: bad" &&
      error.name === "EngineError",
  );
});

test("parseAnswer passes malformed text on as a SyntaxError", () => {
  assert.throws(() => parseAnswer("not json"), SyntaxError);
});

test("wrapEngine passes every argument through in order", () => {
  const calls = [];
  const module = {
    searchUrl: (...args) => (calls.push(["searchUrl", args]), '{"url":"https://x"}'),
    groupNames: (...args) => (calls.push(["groupNames", args]), '["ketone"]'),
    readSpectra: (...args) => (calls.push(["readSpectra", args]), "[]"),
    suggest: (...args) => (calls.push(["suggest", args]), "[]"),
    rank: (...args) => (calls.push(["rank", args]), '{"candidates":[]}'),
  };
  const engine = wrapEngine(module);
  const spectra = { carbon: { format: "pasted", text: "1" }, dept135: null, dept90: null };
  assert.equal(engine.searchUrl("F"), "https://x");
  assert.deepEqual(engine.groupNames(), ["ketone"]);
  assert.deepEqual(engine.readSpectra(spectra), []);
  assert.deepEqual(engine.suggest(spectra, false), []);
  assert.deepEqual(
    engine.rank({
      formula: "F",
      pubchemCsv: "csv",
      spectra,
      present: ["p"],
      absent: ["a"],
      useDept: true,
    }),
    { candidates: [] },
  );
  assert.deepEqual(calls, [
    ["searchUrl", ["F"]],
    ["groupNames", []],
    ["readSpectra", [spectra]],
    ["suggest", [spectra, false]],
    ["rank", ["F", "csv", spectra, ["p"], ["a"], true]],
  ]);
});

test("wrapEngine throws the engine's errors", () => {
  const failing = '{"error":{"code":"InvalidArgument","message":"no"}}';
  const engine = wrapEngine({ searchUrl: () => failing });
  assert.throws(() => engine.searchUrl("F"), { name: "EngineError", code: "InvalidArgument" });
});
