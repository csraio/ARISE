// Fetching a PubChem formula search. The engine builds the address and reads the answer;
// only this module touches the network. `fetchFn` is the browser's fetch, or a stand-in.

/** PubChem could not be reached or refused, with a message for the chemist. */
export class PubChemError extends Error {
  constructor(message) {
    super(message);
    this.name = "PubChemError";
  }
}

/** PubChem's CSV answer for `url`. */
export async function fetchSearch(url, fetchFn = globalThis.fetch) {
  let response;
  try {
    response = await fetchFn(url);
  } catch (error) {
    throw new PubChemError(`could not reach PubChem (${error.message}); check the connection, or load a saved PubChem CSV`);
  }
  if (!response.ok) {
    const status = `${response.status} ${response.statusText}`.trim();
    throw new PubChemError(`PubChem answered ${status}; try again later, or load a saved PubChem CSV`);
  }
  return response.text();
}

/** The number of data rows in a PubChem CSV (lines after the header, blank lines aside). */
export function countRows(csv) {
  const lines = csv.split(/\r?\n/).filter((line) => line !== "");
  return Math.max(0, lines.length - 1);
}
