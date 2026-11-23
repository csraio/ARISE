// public/js/peaks.js: the three ways to enter a spectrum (paste, upload, table) become the
// {format, text} the engine reads. Numbers are checked by the engine, not here.
import assert from "node:assert/strict";
import { test } from "node:test";

import { PeakEntryError, spectrumFrom, tableToCsv } from "../js/peaks.js";

const row = (shift, intensity = "") => ({ shift, intensity });

test("a 13C table without intensities is a shift-only CSV", () => {
  assert.equal(tableToCsv([row("170.5"), row("29")], "carbon"), "shift_ppm\n170.5\n29\n");
});

test("a 13C table with every intensity keeps them", () => {
  assert.equal(
    tableToCsv([row("170.5", "1200"), row("29", "-3")], "carbon"),
    "shift_ppm,intensity\n170.5,1200\n29,-3\n",
  );
});

test("a 13C table with some intensities missing is refused, naming the row", () => {
  assert.throws(
    () => tableToCsv([row("170.5", "1200"), row("29")], "carbon"),
    (error) => error instanceof PeakEntryError && /row 2/.test(error.message),
  );
});

test("a DEPT table needs every intensity", () => {
  assert.equal(
    tableToCsv([row("29", "-850"), row("128", "100")], "dept"),
    "shift_ppm,intensity\n29,-850\n128,100\n",
  );
  assert.throws(() => tableToCsv([row("29", "-850"), row("128")], "dept"), /row 2/);
});

test("empty rows are skipped and cells are trimmed", () => {
  assert.equal(
    tableToCsv([row(" 170.5 "), row(""), row("  ", " "), row("29")], "carbon"),
    "shift_ppm\n170.5\n29\n",
  );
});

test("an intensity without a shift is refused", () => {
  assert.throws(() => tableToCsv([row("170.5"), row("", "40")], "carbon"), /row 2/);
});

test("a comma in a cell is refused (decimals take a point)", () => {
  assert.throws(() => tableToCsv([row("29,5")], "carbon"), /row 1.*point/);
  assert.throws(() => tableToCsv([row("29", "1,000")], "dept"), /row 1/);
});

test("a table without peaks is a header only: the engine decides if that is allowed", () => {
  assert.equal(tableToCsv([], "dept"), "shift_ppm,intensity\n");
  assert.equal(tableToCsv([row("")], "carbon"), "shift_ppm\n");
});

test("spectrumFrom: paste, upload and table", () => {
  assert.deepEqual(spectrumFrom({ mode: "paste", text: "29\n" }, "carbon"), {
    format: "pasted",
    text: "29\n",
  });
  assert.deepEqual(spectrumFrom({ mode: "upload", fileText: "shift_ppm\n29\n" }, "carbon"), {
    format: "csv",
    text: "shift_ppm\n29\n",
  });
  assert.deepEqual(spectrumFrom({ mode: "table", rows: [row("29", "-1")] }, "dept"), {
    format: "csv",
    text: "shift_ppm,intensity\n29,-1\n",
  });
});

test("spectrumFrom: a DEPT spectrum may be not recorded; 13C may not", () => {
  assert.equal(spectrumFrom({ mode: "none" }, "dept"), null);
  assert.throws(() => spectrumFrom({ mode: "none" }, "carbon"), PeakEntryError);
});

test("spectrumFrom: an upload without a file is refused", () => {
  assert.throws(() => spectrumFrom({ mode: "upload", fileText: null }, "carbon"), /file/);
});

test("spectrumFrom: an unknown mode is a programming error", () => {
  assert.throws(() => spectrumFrom({ mode: "fax" }, "carbon"), /unknown/);
});
