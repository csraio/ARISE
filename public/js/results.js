// Summaries of the engine's ranking, and the ranking as CSV for a spreadsheet.

/** How many candidates each tier holds. */
export function tierCounts(ranking) {
  const counts = { active: 0, "back-burner": 0, eliminated: 0 };
  for (const candidate of ranking.candidates) counts[candidate.tier] += 1;
  return counts;
}

/** How many PubChem records cleaning removed, by reason. */
export function removalCounts(ranking) {
  const counts = {};
  for (const removal of ranking.removed) counts[removal.reason] = (counts[removal.reason] ?? 0) + 1;
  return counts;
}

function csvField(value) {
  const text = String(value);
  return /[",\n\r]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

/** The ranking as CSV, one row per candidate in rank order; reasons and notes joined by "; ". */
export function rankingCsv(ranking) {
  const header = "rank,cid,smiles,tier,assumed_overlaps,environments_fewest,environments_most,reasons,notes\n";
  return (
    header +
    ranking.candidates
      .map((c, i) =>
        [
          i + 1,
          c.cid,
          c.smiles,
          c.tier,
          c.assumed_overlaps,
          c.environments_fewest,
          c.environments_most,
          c.reasons.join("; "),
          c.notes.join("; "),
        ]
          .map(csvField)
          .join(",") + "\n",
      )
      .join("")
  );
}

/** The PubChem page of a compound. */
export function pubchemLink(cid) {
  return `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;
}
