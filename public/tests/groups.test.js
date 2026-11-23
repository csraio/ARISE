// public/js/groups.js: which groups the page offers, in which order, and what the chemist's
// ticks mean for the engine.
import assert from "node:assert/strict";
import { test } from "node:test";

import { confirmedFrom, groupChoices } from "../js/groups.js";

const all = ["acid", "ketone", "ester", "ether"];
const suggestions = [
  { group: "ester", evidence: [{ spectrum: "13C", peak: 0 }] },
  { group: "ketone", evidence: [{ spectrum: "13C", peak: 1 }, { spectrum: "1H", peak: 2 }] },
];

test("suggested groups come first, in suggestion order, with their evidence", () => {
  assert.deepEqual(groupChoices(all, suggestions), [
    { group: "ester", suggested: true, fromProton: false, evidence: [{ spectrum: "13C", peak: 0 }] },
    {
      group: "ketone",
      suggested: true,
      fromProton: true,
      evidence: [{ spectrum: "13C", peak: 1 }, { spectrum: "1H", peak: 2 }],
    },
    { group: "acid", suggested: false, fromProton: false, evidence: [] },
    { group: "ether", suggested: false, fromProton: false, evidence: [] },
  ]);
});

test("without suggestions every group is offered in engine order", () => {
  assert.deepEqual(
    groupChoices(all, []).map((c) => c.group),
    all,
  );
});

test("a suggestion for a group the engine does not list is a programming error", () => {
  assert.throws(() => groupChoices(all, [{ group: "nitrile", evidence: [] }]), /nitrile/);
});

test("ticks become present, absent and likely lists in engine order; unsure is left out", () => {
  const ticks = { ether: "absent", acid: "present", ketone: "likely", ester: "present" };
  assert.deepEqual(confirmedFrom(all, ticks), {
    present: ["acid", "ester"],
    absent: ["ether"],
    likely: ["ketone"],
  });
  assert.deepEqual(confirmedFrom(all, { ketone: "unsure" }), { present: [], absent: [], likely: [] });
  assert.deepEqual(confirmedFrom(all, {}), { present: [], absent: [], likely: [] });
});

test("an unknown tick value is a programming error", () => {
  assert.throws(() => confirmedFrom(all, { acid: "maybe" }), /maybe/);
});
