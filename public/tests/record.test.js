// public/js/record.js: the downloadable record of a search (formula, query date, the list
// PubChem returned, the spectra and the chemist's choices) and its re-import.
import assert from "node:assert/strict";
import { test } from "node:test";

import { RecordError, makeRecord, readRecord, recordFileName, recordJson, sameRanking } from "../js/record.js";

const inputs = {
  formula: "C3H6O",
  queryDate: "2026-11-22T18:30:00.000Z",
  pubchemUrl: "https://pubchem.example/x",
  pubchemCsv: '"CID","ConnectivitySMILES"\n',
  spectra: {
    carbon: { format: "pasted", text: "29\n" },
    dept135: { format: "csv", text: "shift_ppm,intensity\n29,-1\n" },
    dept90: null,
  },
  present: ["ketone"],
  absent: [],
  useDept: true,
  ranking: { formula: "C3H6O", rows: 0, removed: [], candidates: [] },
};

test("a record reads back as written", () => {
  const record = makeRecord(inputs);
  assert.equal(record.format, "arise-record");
  assert.equal(record.version, 1);
  assert.deepEqual(readRecord(recordJson(record)), record);
});

test("a record without a ranking is valid", () => {
  const record = makeRecord({ ...inputs, ranking: null });
  assert.deepEqual(readRecord(recordJson(record)), record);
});

test("the file name carries the formula and the query day", () => {
  assert.equal(recordFileName("C3H6O", "2026-11-22T18:30:00.000Z"), "arise-C3H6O-2026-11-22.json");
});

function rejects(text, pattern) {
  assert.throws(() => readRecord(text), (error) => error instanceof RecordError && pattern.test(error.message));
}

test("text that is not a record is refused with a reason", () => {
  rejects("not json", /not an ARISE record/);
  rejects("[]", /not an ARISE record/);
  rejects(JSON.stringify({ format: "something else", version: 1 }), /not an ARISE record/);
  rejects(JSON.stringify({ ...makeRecord(inputs), version: 2 }), /version 2/);
});

test("a record with a missing or wrong field is refused, naming the field", () => {
  const record = makeRecord(inputs);
  for (const [field, value] of [
    ["formula", 7],
    ["queryDate", null],
    ["pubchemUrl", undefined],
    ["pubchemCsv", 5],
    ["useDept", "yes"],
    ["present", "ketone"],
    ["absent", [3]],
    ["spectra", null],
  ]) {
    const broken = { ...record, [field]: value };
    rejects(JSON.stringify(broken), new RegExp(field));
  }
  rejects(JSON.stringify({ ...record, spectra: { ...record.spectra, carbon: null } }), /carbon/);
  rejects(
    JSON.stringify({ ...record, spectra: { ...record.spectra, dept90: { format: "xml", text: "" } } }),
    /dept90/,
  );
});

test("sameRanking compares whole rankings", () => {
  const a = { formula: "F", rows: 1, removed: [], candidates: [{ cid: 1, tier: "active" }] };
  assert.ok(sameRanking(a, structuredClone(a)));
  assert.ok(!sameRanking(a, { ...a, rows: 2 }));
  assert.ok(!sameRanking(a, { ...a, candidates: [{ cid: 1, tier: "eliminated" }] }));
});
