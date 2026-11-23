// public/js/results.js: summaries of a ranking and its CSV for spreadsheets.
import assert from "node:assert/strict";
import { test } from "node:test";

import { pubchemLink, rankingCsv, removalCounts, tierCounts } from "../js/results.js";

const ranking = {
  formula: "F",
  rows: 5,
  removed: [
    { cid: 9, smiles: "x", reason: "mixture" },
    { cid: 8, smiles: "y", reason: "mixture" },
    { cid: 7, smiles: "z", reason: "charged" },
  ],
  candidates: [
    {
      cid: 1,
      smiles: "first-structure",
      tier: "active",
      reasons: [],
      notes: ["assumes 1 overlapping peak"],
      assumed_overlaps: 1,
      environments_fewest: 3,
      environments_most: 4,
    },
    {
      cid: 2,
      smiles: "second-structure",
      tier: "eliminated",
      reasons: ['lacks "ketone"', "second, reason"],
      notes: [],
      assumed_overlaps: 0,
      environments_fewest: 0,
      environments_most: 0,
    },
  ],
};

test("tierCounts counts every tier, zero included", () => {
  assert.deepEqual(tierCounts(ranking), { active: 1, "back-burner": 0, eliminated: 1 });
});

test("removalCounts counts each reason", () => {
  assert.deepEqual(removalCounts(ranking), { mixture: 2, charged: 1 });
  assert.deepEqual(removalCounts({ ...ranking, removed: [] }), {});
});

test("rankingCsv has one row per candidate in rank order, quoted where needed", () => {
  assert.equal(
    rankingCsv(ranking),
    "rank,cid,smiles,tier,assumed_overlaps,environments_fewest,environments_most,reasons,notes\n" +
      "1,1,first-structure,active,1,3,4,,assumes 1 overlapping peak\n" +
      '2,2,second-structure,eliminated,0,0,0,"lacks ""ketone""; second, reason",\n',
  );
});

test("pubchemLink points at the compound page", () => {
  assert.equal(pubchemLink(2244), "https://pubchem.ncbi.nlm.nih.gov/compound/2244");
});
