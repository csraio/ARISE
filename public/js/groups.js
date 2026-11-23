// The functional groups the page offers for the chemist to confirm. Suggestions never
// eliminate; only a group ticked present or absent is passed to the engine.

/**
 * Every group to show: those the spectra suggest first (in the engine's suggestion order,
 * with their evidence), then the rest in the engine's order. `fromProton` marks groups with
 * 1H evidence: the chemist may confirm those as "likely" (soft, the chemist's decision
 * 2026-11-22).
 */
export function groupChoices(allNames, suggestions) {
  const known = new Set(allNames);
  const suggested = new Map();
  for (const suggestion of suggestions) {
    if (!known.has(suggestion.group)) throw new Error(`the engine suggested an unknown group "${suggestion.group}"`);
    suggested.set(suggestion.group, suggestion.evidence);
  }
  return [
    ...[...suggested].map(([group, evidence]) => ({
      group,
      suggested: true,
      fromProton: evidence.some((e) => e.spectrum === "1H"),
      evidence,
    })),
    ...allNames
      .filter((group) => !suggested.has(group))
      .map((group) => ({ group, suggested: false, fromProton: false, evidence: [] })),
  ];
}

const kTicks = new Set(["present", "absent", "likely", "unsure"]);

/**
 * The engine's confirmed groups from the ticks {group: "present" | "absent" | "likely" |
 * "unsure"}: present and absent eliminate, likely back-burners.
 */
export function confirmedFrom(allNames, ticks) {
  for (const value of Object.values(ticks)) {
    if (!kTicks.has(value)) throw new Error(`unknown tick "${value}"`);
  }
  return {
    present: allNames.filter((group) => ticks[group] === "present"),
    absent: allNames.filter((group) => ticks[group] === "absent"),
    likely: allNames.filter((group) => ticks[group] === "likely"),
  };
}
