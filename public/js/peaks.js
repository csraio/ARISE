// The three ways to enter a spectrum (the chemist's choice, 2026-11-22): paste a list, upload
// a CSV file, or type into a table. Each becomes the {format, text} the engine reads. The
// numbers themselves are checked by the engine, which names the bad line.

/** A table or entry the page cannot turn into text, with a message for the chemist. */
export class PeakEntryError extends Error {
  constructor(message) {
    super(message);
    this.name = "PeakEntryError";
  }
}

/**
 * The CSV text of table rows [{shift, intensity}] in ARISE's peak-list format. Empty rows are
 * skipped. `kind` "carbon": intensities are optional, but all or none. `kind` "dept": every
 * row needs one (its sign is the phase).
 */
export function tableToCsv(rows, kind) {
  const filled = [];
  rows.forEach((row, index) => {
    const shift = String(row.shift ?? "").trim();
    const intensity = String(row.intensity ?? "").trim();
    if (shift === "" && intensity === "") return;
    const number = index + 1;
    if (shift === "") throw new PeakEntryError(`row ${number} has an intensity but no shift`);
    if (shift.includes(",") || intensity.includes(",")) {
      throw new PeakEntryError(`row ${number}: use a point, not a comma, for decimals`);
    }
    filled.push({ number, shift, intensity });
  });
  const withIntensity = filled.filter((row) => row.intensity !== "");
  const allIntensities = kind === "dept" || (filled.length > 0 && withIntensity.length === filled.length);
  if (allIntensities) {
    const missing = filled.find((row) => row.intensity === "");
    if (missing !== undefined) {
      throw new PeakEntryError(
        kind === "dept"
          ? `row ${missing.number} has no intensity (DEPT needs it: its sign is the phase)`
          : `row ${missing.number} has no intensity; give every row an intensity or none`,
      );
    }
    return "shift_ppm,intensity\n" + filled.map((r) => `${r.shift},${r.intensity}\n`).join("");
  }
  if (withIntensity.length > 0) {
    const missing = filled.find((row) => row.intensity === "");
    throw new PeakEntryError(`row ${missing.number} has no intensity; give every row an intensity or none`);
  }
  return "shift_ppm\n" + filled.map((r) => `${r.shift}\n`).join("");
}

/**
 * The engine's {format, text} for one spectrum entry, or null for a DEPT spectrum that was
 * not recorded. `entry.mode` is "paste" (entry.text), "upload" (entry.fileText, the file's
 * text or null), "table" (entry.rows) or "none" (DEPT only).
 */
export function spectrumFrom(entry, kind) {
  switch (entry.mode) {
    case "paste":
      return { format: "pasted", text: entry.text ?? "" };
    case "upload":
      if (entry.fileText === null || entry.fileText === undefined) {
        throw new PeakEntryError("choose a CSV file to upload");
      }
      return { format: "csv", text: entry.fileText };
    case "table":
      return { format: "csv", text: tableToCsv(entry.rows ?? [], kind) };
    case "none":
      if (kind === "carbon") throw new PeakEntryError("a 13C spectrum is required");
      return null;
    default:
      throw new Error(`unknown entry mode "${entry.mode}"`);
  }
}
