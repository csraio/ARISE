// The ARISE engine as the page uses it. The engine (apps/wasm) answers every call with JSON
// text; this module parses it and turns an {"error": {code, message}} answer into a thrown
// EngineError, so the page handles failures in one place.
//
// The engine module is passed in (`createAriseEngine` from arise_engine.js), so the page and
// the Node tests load it from their own locations.

/** A failure reported by the engine, with its code ("ParseError", "InvalidArgument", ...). */
export class EngineError extends Error {
  constructor(code, message) {
    super(message);
    this.name = "EngineError";
    this.code = code;
  }
}

/** Parses one engine answer; throws EngineError for an error answer. */
export function parseAnswer(text) {
  const answer = JSON.parse(text);
  if (answer !== null && typeof answer === "object" && "error" in answer) {
    throw new EngineError(answer.error.code, answer.error.message);
  }
  return answer;
}

/**
 * Wraps an instantiated engine module (the object createAriseEngine resolves to).
 *
 * Spectra are {carbon, dept135, dept90}; each spectrum is {format: "csv" | "pasted", text},
 * and dept135/dept90 are null when not recorded.
 */
export function wrapEngine(module) {
  return {
    /** The PubChem address listing every compound with `formula`. */
    searchUrl(formula) {
      return parseAnswer(module.searchUrl(formula)).url;
    },
    /** Every functional group's name, in the engine's order. */
    groupNames() {
      return parseAnswer(module.groupNames());
    },
    /** The 13C peaks with their possible DEPT types: [{shift, types}]. */
    readSpectra(spectra) {
      return parseAnswer(module.readSpectra(spectra));
    },
    /** Group suggestions with their evidence: [{group, evidence: [...]}]. */
    suggest(spectra, useDept) {
      return parseAnswer(module.suggest(spectra, useDept));
    },
    /** The cleaned PubChem list, ranked: {formula, rows, removed, candidates}. */
    rank({ formula, pubchemCsv, spectra, present, absent, useDept }) {
      return parseAnswer(module.rank(formula, pubchemCsv, spectra, present, absent, useDept));
    },
  };
}

/** Instantiates the engine module and wraps it. */
export async function loadEngine(createAriseEngine) {
  return wrapEngine(await createAriseEngine());
}
