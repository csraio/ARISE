// The record of one search, to download and re-import: the formula, when PubChem was asked,
// the list it returned, the spectra as entered, and the chemist's choices. PubChem changes
// over time, so the list itself is kept: re-importing gives the same input, and so the
// same ranking (golden rule 4).

export const kRecordFormat = "arise-record";
export const kRecordVersion = 2;  // 2: adds the 1H spectrum and the likely groups

/** A file that is not a usable record, with the reason. */
export class RecordError extends Error {
  constructor(message) {
    super(message);
    this.name = "RecordError";
  }
}

/** A record from the page's state; `ranking` is the engine's answer, or null. */
export function makeRecord({ formula, queryDate, pubchemUrl, pubchemCsv, spectra, present, absent, likely, useDept, ranking }) {
  return {
    format: kRecordFormat,
    version: kRecordVersion,
    formula,
    queryDate,
    pubchemUrl,
    pubchemCsv,
    spectra,
    present,
    absent,
    likely: likely ?? [],
    useDept,
    ranking: ranking ?? null,
  };
}

/** The record as the text of a .json file. */
export function recordJson(record) {
  return JSON.stringify(record, null, 1) + "\n";
}

/** The download name: arise-<formula>-<query day>.json. */
export function recordFileName(formula, queryDate) {
  return `arise-${formula}-${queryDate.slice(0, 10)}.json`;
}

function requireType(record, field, type) {
  if (typeof record[field] !== type) throw new RecordError(`the record's ${field} is missing or not ${type === "string" ? "text" : type}`);
}

function requireNames(record, field) {
  if (!Array.isArray(record[field]) || !record[field].every((name) => typeof name === "string")) {
    throw new RecordError(`the record's ${field} must be a list of group names`);
  }
}

function requireSpectrum(spectra, field, optional) {
  const spectrum = spectra[field];
  if (optional && spectrum === null) return;
  if (
    spectrum === null ||
    typeof spectrum !== "object" ||
    !["csv", "pasted"].includes(spectrum.format) ||
    typeof spectrum.text !== "string"
  ) {
    throw new RecordError(`the record's spectra.${field} is missing or malformed`);
  }
}

/** Reads and checks a record file's text; throws RecordError with the reason. */
export function readRecord(text) {
  let record;
  try {
    record = JSON.parse(text);
  } catch {
    throw new RecordError("this file is not an ARISE record (it is not JSON)");
  }
  if (record === null || typeof record !== "object" || Array.isArray(record) || record.format !== kRecordFormat) {
    throw new RecordError("this file is not an ARISE record");
  }
  if (record.version === 1) {
    // Version 1 had no 1H spectrum and no likely groups.
    record = { ...record, version: kRecordVersion, likely: [], spectra: { ...record.spectra, proton: null } };
  }
  if (record.version !== kRecordVersion) {
    throw new RecordError(`this record is version ${record.version}; this page reads versions 1 and ${kRecordVersion}`);
  }
  for (const field of ["formula", "queryDate", "pubchemUrl", "pubchemCsv"]) requireType(record, field, "string");
  requireType(record, "useDept", "boolean");
  requireNames(record, "present");
  requireNames(record, "absent");
  requireNames(record, "likely");
  if (record.spectra === null || typeof record.spectra !== "object") {
    throw new RecordError("the record's spectra are missing");
  }
  requireSpectrum(record.spectra, "carbon", false);
  requireSpectrum(record.spectra, "dept135", true);
  requireSpectrum(record.spectra, "dept90", true);
  requireSpectrum(record.spectra, "proton", true);
  if (record.ranking !== null && (typeof record.ranking !== "object" || !Array.isArray(record.ranking.candidates))) {
    throw new RecordError("the record's ranking is malformed");
  }
  return record;
}

/** Whether two engine rankings are identical (both come from the engine's fixed-order JSON). */
export function sameRanking(a, b) {
  return JSON.stringify(a) === JSON.stringify(b);
}
