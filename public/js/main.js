// The page: connects the controls to the engine and to the tested modules (engine.js,
// peaks.js, pubchem.js, groups.js, results.js, record.js). Everything chemical happens in the
// engine; this file only moves values between the page and those modules.

import createAriseEngine from "../engine/arise_engine.js";
import { EngineError, loadEngine } from "./engine.js";
import { confirmedFrom, groupChoices } from "./groups.js";
import { PeakEntryError, spectrumFrom } from "./peaks.js";
import { PubChemError, countRows, fetchSearch } from "./pubchem.js";
import { RecordError, makeRecord, readRecord, recordFileName, recordJson, sameRanking } from "./record.js";
import { pubchemLink, rankingCsv, removalCounts, tierCounts } from "./results.js";

const kPageSize = 100;
const kTableRows = 6;

const $ = (id) => document.getElementById(id);

const state = {
  engine: null,
  groupNames: [],
  search: null, // {formula, queryDate, pubchemUrl, pubchemCsv}
  suggestions: [],
  ticks: {},
  ranking: null,
  shown: 0,
};

// ---- messages ----------------------------------------------------------------------

function showError(id, error) {
  const box = $(id);
  if (error === null) {
    box.hidden = true;
    box.textContent = "";
    return;
  }
  const known = [EngineError, PeakEntryError, PubChemError, RecordError].some((type) => error instanceof type);
  box.textContent = known ? error.message : `Unexpected problem: ${error.message}`;
  box.hidden = false;
  if (!known) console.error(error);
}

function setStatus(id, text) {
  $(id).textContent = text;
}

// ---- 1. candidates -----------------------------------------------------------------

function useSearch(search, source) {
  state.search = search;
  state.ranking = null;
  $("ranking").hidden = true;
  setStatus("ranking-note", "");
  const rows = countRows(search.pubchemCsv);
  const day = new Date(search.queryDate).toLocaleString();
  setStatus("candidates-status", `${search.formula}: ${rows} PubChem records (${source}, ${day}).`);
}

async function searchPubChem(event) {
  event.preventDefault();
  showError("candidates-error", null);
  const formula = $("formula").value.trim();
  const button = $("search-button");
  button.disabled = true;
  setStatus("candidates-status", "Asking PubChem…");
  try {
    const pubchemUrl = state.engine.searchUrl(formula);
    const pubchemCsv = await fetchSearch(pubchemUrl);
    useSearch({ formula, queryDate: new Date().toISOString(), pubchemUrl, pubchemCsv }, "PubChem");
  } catch (error) {
    setStatus("candidates-status", "");
    showError("candidates-error", error);
  } finally {
    button.disabled = false;
  }
}

async function loadPubChemFile(event) {
  showError("candidates-error", null);
  const file = event.target.files[0];
  if (file === undefined) return;
  try {
    const formula = $("formula").value.trim();
    if (formula === "") throw new PeakEntryError("type the molecular formula the file was searched for first");
    // Checks the formula; the address records where such a list comes from.
    const pubchemUrl = state.engine.searchUrl(formula);
    const pubchemCsv = await file.text();
    // The query date of a saved list is not known; the file's date is the closest.
    const queryDate = new Date(file.lastModified).toISOString();
    useSearch({ formula, queryDate, pubchemUrl, pubchemCsv }, `file ${file.name}`);
  } catch (error) {
    showError("candidates-error", error);
  } finally {
    event.target.value = "";
  }
}

// ---- 2. spectra --------------------------------------------------------------------

const spectra = {}; // name -> {root, kind, fileText}

function addTableRow(tbody, shift = "", intensity = "") {
  const tr = document.createElement("tr");
  for (const [field, value] of [["shift", shift], ["intensity", intensity]]) {
    const td = document.createElement("td");
    const input = document.createElement("input");
    input.inputMode = "decimal";
    input.dataset.field = field;
    input.value = value;
    input.setAttribute("aria-label", field === "shift" ? "Shift (ppm)" : "Intensity");
    td.append(input);
    tr.append(td);
  }
  tbody.append(tr);
}

function showMode(spectrum) {
  const mode = spectrum.root.querySelector(".mode-select").value;
  for (const entry of spectrum.root.querySelectorAll(".entry")) entry.hidden = entry.dataset.mode !== mode;
}

function buildSpectrumPanels() {
  const template = $("spectrum-template");
  for (const holder of document.querySelectorAll("[data-spectrum]")) {
    const root = template.content.firstElementChild.cloneNode(true);
    const name = holder.dataset.spectrum;
    const kind = holder.dataset.kind;
    root.querySelector("legend").textContent = holder.dataset.title;
    const select = root.querySelector(".mode-select");
    if (kind === "proton") {
      // The table's second column is the integral, as a number of protons.
      root.querySelectorAll(".entry-table th")[1].textContent = "Protons (integral)";
    }
    if (kind === "carbon") {
      select.querySelector('option[value="none"]').remove();
      select.value = "paste";
    }
    const tbody = root.querySelector(".entry-table tbody");
    for (let i = 0; i < kTableRows; i++) addTableRow(tbody);
    root.querySelector(".add-row").addEventListener("click", () => addTableRow(tbody));
    const spectrum = { root, kind, fileText: null };
    root.querySelector('.entry[data-mode="upload"] input').addEventListener("change", async (event) => {
      const file = event.target.files[0];
      spectrum.fileText = file === undefined ? null : await file.text();
      root.querySelector(".file-status").textContent = file === undefined ? "" : `Loaded ${file.name}`;
    });
    select.addEventListener("change", () => showMode(spectrum));
    holder.append(root);
    spectra[name] = spectrum;
    showMode(spectrum);
  }
}

function entryOf(spectrum) {
  const root = spectrum.root;
  return {
    mode: root.querySelector(".mode-select").value,
    text: root.querySelector("textarea").value,
    fileText: spectrum.fileText,
    rows: [...root.querySelectorAll(".entry-table tbody tr")].map((tr) => ({
      shift: tr.querySelector('[data-field="shift"]').value,
      intensity: tr.querySelector('[data-field="intensity"]').value,
    })),
  };
}

/** The engine's spectra from the three panels; throws PeakEntryError naming the panel. */
function currentSpectra() {
  const result = {};
  for (const [name, spectrum] of Object.entries(spectra)) {
    try {
      result[name] = spectrumFrom(entryOf(spectrum), spectrum.kind);
    } catch (error) {
      if (error instanceof PeakEntryError) {
        throw new PeakEntryError(`${spectrum.root.querySelector("legend").textContent}: ${error.message}`);
      }
      throw error;
    }
  }
  return result;
}

/** Puts a saved spectrum back into its panel (pasted text, or CSV as an upload). */
function restoreSpectrum(name, saved) {
  const spectrum = spectra[name];
  const select = spectrum.root.querySelector(".mode-select");
  spectrum.fileText = null;
  spectrum.root.querySelector(".file-status").textContent = "";
  if (saved === null) {
    select.value = "none";
  } else if (saved.format === "pasted") {
    select.value = "paste";
    spectrum.root.querySelector("textarea").value = saved.text;
  } else {
    select.value = "upload";
    spectrum.fileText = saved.text;
    spectrum.root.querySelector(".file-status").textContent = "Loaded from the saved results";
  }
  showMode(spectrum);
}

function showPeaks(peaks) {
  const tbody = $("peaks-table").querySelector("tbody");
  tbody.replaceChildren(
    ...peaks.map((peak, i) => {
      const tr = document.createElement("tr");
      for (const [text, className] of [
        [String(i + 1), "number"],
        [String(peak.shift), "number"],
        [peak.types.join(", "), ""],
      ]) {
        const td = document.createElement("td");
        td.textContent = text;
        if (className !== "") td.className = className;
        tr.append(td);
      }
      return tr;
    }),
  );
  $("peaks").hidden = false;
}

function showProtons(protons) {
  const box = $("protons");
  if (protons === null) {
    box.hidden = true;
    return;
  }
  const tbody = $("protons-table").querySelector("tbody");
  tbody.replaceChildren(
    ...protons.map((peak, i) => {
      const tr = document.createElement("tr");
      for (const text of [String(i + 1), String(peak.shift), peak.protons === null ? "–" : `${peak.protons}H`]) {
        const td = document.createElement("td");
        td.textContent = text;
        td.className = "number";
        tr.append(td);
      }
      return tr;
    }),
  );
  box.hidden = false;
}

function readSpectra() {
  showError("spectra-error", null);
  try {
    const current = currentSpectra();
    showPeaks(state.engine.readSpectra(current));
    showProtons(state.engine.readProtons(current));
    suggestGroups(current);
  } catch (error) {
    $("peaks").hidden = true;
    $("protons").hidden = true;
    showError("spectra-error", error);
  }
}

// ---- 3. groups ---------------------------------------------------------------------

const percent = (rate) => `${Math.round(rate * 100)}%`;

function evidenceText(e) {
  const where = `${e.low}–${e.high} ppm`;
  if (e.spectrum === "1H") {
    const protons = e.protons === null ? "" : ` (${e.protons}H)`;
    return (
      `1H peak ${e.peak + 1} at ${e.shift} ppm${protons}: H on the ${e.hallmark} (${e.h_count}) window ${where}; ` +
      `${percent(e.precision)} of nmrshiftdb2 molecules with such a peak have the group ` +
      `(${percent(e.base_rate)} of all molecules do)`
    );
  }
  const basis = e.source === "lab" ? "your lab's window" : `${e.molecules} molecules, ${e.source}`;
  return `13C peak ${e.peak + 1} at ${e.shift} ppm: ${e.hallmark} window ${where} (${basis})`;
}

function groupItem(choice) {
  const item = document.createElement("div");
  item.className = choice.suggested ? "group suggested" : "group";
  const head = document.createElement("div");
  head.className = "group-head";
  const name = document.createElement("span");
  name.className = "group-name";
  name.textContent = choice.group;
  const choices = document.createElement("fieldset");
  choices.className = "choices";
  choices.setAttribute("aria-label", choice.group);
  // "likely" only for groups with 1H evidence (candidates without them are back-burnered), or
  // when it is already chosen (a reopened record), so the choice in effect is always shown.
  const offerLikely = choice.fromProton || state.ticks[choice.group] === "likely";
  const values = offerLikely ? ["unsure", "likely", "present", "absent"] : ["unsure", "present", "absent"];
  for (const value of values) {
    const label = document.createElement("label");
    const radio = document.createElement("input");
    radio.type = "radio";
    radio.name = `group-${choice.group}`;
    radio.value = value;
    radio.checked = (state.ticks[choice.group] ?? "unsure") === value;
    radio.addEventListener("change", () => {
      state.ticks[choice.group] = value;
    });
    label.append(radio, ` ${value}`);
    choices.append(label);
  }
  head.append(name, choices);
  item.append(head);
  if (choice.evidence.length > 0) {
    const list = document.createElement("ul");
    list.className = "evidence";
    for (const e of choice.evidence) {
      const li = document.createElement("li");
      li.textContent = evidenceText(e);
      list.append(li);
    }
    item.append(list);
  }
  return item;
}

function showGroups() {
  const choices = groupChoices(state.groupNames, state.suggestions);
  const suggested = choices.filter((c) => c.suggested);
  const others = choices.filter((c) => !c.suggested);
  const box = $("groups");
  const intro = document.createElement("p");
  intro.className = "hint";
  intro.textContent =
    suggested.length === 0
      ? "The spectra suggest no group. You can still mark any group below."
      : `The spectra suggest ${suggested.length} group${suggested.length === 1 ? "" : "s"}:`;
  const more = document.createElement("details");
  more.className = "more-groups";
  const summary = document.createElement("summary");
  summary.textContent = `Other groups (${others.length})`;
  more.append(summary, ...others.map(groupItem));
  more.open = others.some((c) => (state.ticks[c.group] ?? "unsure") !== "unsure");
  box.replaceChildren(intro, ...suggested.map(groupItem), more);
}

function suggestGroups(current) {
  showError("groups-error", null);
  try {
    state.suggestions = state.engine.suggest(current, $("use-dept").checked);
    showGroups();
  } catch (error) {
    showError("groups-error", error);
  }
}

// ---- 4. ranking --------------------------------------------------------------------

function showSummary(ranking) {
  const tiers = tierCounts(ranking);
  const removed = removalCounts(ranking);
  const items = [
    ["PubChem records", ranking.rows],
    ["Distinct candidates", ranking.candidates.length],
    ["Active", tiers.active],
    ["Back-burner", tiers["back-burner"]],
    ["Eliminated", tiers.eliminated],
    ...Object.entries(removed).map(([reason, n]) => [`Removed: ${reason}`, n]),
  ];
  $("ranking-summary").replaceChildren(
    ...items.map(([label, value]) => {
      const div = document.createElement("div");
      const dt = document.createElement("dt");
      dt.textContent = label;
      const dd = document.createElement("dd");
      dd.textContent = String(value);
      div.append(dt, dd);
      return div;
    }),
  );
}

function candidateRow(candidate, rank) {
  const tr = document.createElement("tr");
  const cell = (className = "") => {
    const td = document.createElement("td");
    if (className !== "") td.className = className;
    tr.append(td);
    return td;
  };
  cell("number").textContent = String(rank);
  const link = document.createElement("a");
  link.href = pubchemLink(candidate.cid);
  link.target = "_blank";
  link.rel = "noopener";
  link.textContent = String(candidate.cid);
  cell("number").append(link);
  const tier = cell(`tier tier-${candidate.tier}`);
  tier.textContent = candidate.tier;
  cell().textContent =
    candidate.environments_most === 0
      ? "–"
      : candidate.environments_fewest === candidate.environments_most
        ? String(candidate.environments_most)
        : `${candidate.environments_fewest}–${candidate.environments_most}`;
  const notes = [...candidate.reasons, ...candidate.notes];
  const reasons = cell();
  if (notes.length > 0) {
    const list = document.createElement("ul");
    for (const note of notes) {
      const li = document.createElement("li");
      li.textContent = note;
      list.append(li);
    }
    reasons.append(list);
  }
  cell("smiles").textContent = candidate.smiles;
  return tr;
}

function filteredCandidates() {
  const tier = $("tier-filter").value;
  return state.ranking.candidates
    .map((candidate, i) => ({ candidate, rank: i + 1 }))
    .filter(({ candidate }) => tier === "all" || candidate.tier === tier);
}

function showRows(reset) {
  const tbody = $("ranking-table").querySelector("tbody");
  if (reset) {
    tbody.replaceChildren();
    state.shown = 0;
  }
  const rows = filteredCandidates();
  const next = rows.slice(state.shown, state.shown + kPageSize);
  tbody.append(...next.map(({ candidate, rank }) => candidateRow(candidate, rank)));
  state.shown += next.length;
  const more = $("more-button");
  more.hidden = state.shown >= rows.length;
  more.textContent = `Show more (${rows.length - state.shown} left)`;
}

function showRanking(ranking) {
  state.ranking = ranking;
  showSummary(ranking);
  showRows(true);
  $("ranking").hidden = false;
}

function rankCandidates() {
  showError("ranking-error", null);
  setStatus("ranking-note", "");
  try {
    if (state.search === null) throw new PeakEntryError("search PubChem or load a list first (step 1)");
    const ranking = state.engine.rank({
      formula: state.search.formula,
      pubchemCsv: state.search.pubchemCsv,
      spectra: currentSpectra(),
      ...confirmedFrom(state.groupNames, state.ticks),
      useDept: $("use-dept").checked,
    });
    showRanking(ranking);
  } catch (error) {
    showError("ranking-error", error);
  }
}

// ---- downloads and re-import -------------------------------------------------------

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  document.body.append(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function downloadRecord() {
  showError("ranking-error", null);
  try {
    const record = makeRecord({
      ...state.search,
      spectra: currentSpectra(),
      ...confirmedFrom(state.groupNames, state.ticks),
      useDept: $("use-dept").checked,
      ranking: state.ranking,
    });
    download(recordFileName(record.formula, record.queryDate), recordJson(record), "application/json");
  } catch (error) {
    showError("ranking-error", error);
  }
}

function downloadCsv() {
  const name = recordFileName(state.search.formula, state.search.queryDate).replace(/\.json$/, "-ranking.csv");
  download(name, rankingCsv(state.ranking), "text/csv");
}

async function openRecord(event) {
  showError("candidates-error", null);
  const file = event.target.files[0];
  if (file === undefined) return;
  try {
    const record = readRecord(await file.text());
    $("formula").value = record.formula;
    useSearch(
      { formula: record.formula, queryDate: record.queryDate, pubchemUrl: record.pubchemUrl, pubchemCsv: record.pubchemCsv },
      `saved results ${file.name}`,
    );
    for (const name of Object.keys(spectra)) restoreSpectrum(name, record.spectra[name]);
    $("use-dept").checked = record.useDept;
    state.ticks = {};
    for (const group of record.present) state.ticks[group] = "present";
    for (const group of record.absent) state.ticks[group] = "absent";
    for (const group of record.likely) state.ticks[group] = "likely";
    readSpectra();
    rankCandidates();
    if (record.ranking !== null && state.ranking !== null) {
      setStatus(
        "ranking-note",
        sameRanking(record.ranking, state.ranking)
          ? "Re-ranked: the same result as the saved one."
          : "Re-ranked: the result differs from the saved one (this version of ARISE ranks differently).",
      );
    }
  } catch (error) {
    showError("candidates-error", error);
  } finally {
    event.target.value = "";
  }
}

// ---- start -------------------------------------------------------------------------

async function start() {
  buildSpectrumPanels();
  try {
    state.engine = await loadEngine(createAriseEngine);
    state.groupNames = state.engine.groupNames();
    setStatus("engine-status", "");
  } catch (error) {
    setStatus("engine-status", `The engine did not load: ${error.message}`);
    for (const button of document.querySelectorAll("button")) button.disabled = true;
    return;
  }
  $("search-form").addEventListener("submit", searchPubChem);
  $("pubchem-file").addEventListener("change", loadPubChemFile);
  $("record-file").addEventListener("change", openRecord);
  $("check-button").addEventListener("click", readSpectra);
  $("use-dept").addEventListener("change", () => {
    try {
      suggestGroups(currentSpectra());
    } catch {
      // The spectra are not readable yet; suggestions come with "Read spectra".
    }
  });
  $("rank-button").addEventListener("click", rankCandidates);
  $("tier-filter").addEventListener("change", () => showRows(true));
  $("more-button").addEventListener("click", () => showRows(false));
  $("download-record").addEventListener("click", downloadRecord);
  $("download-csv").addEventListener("click", downloadCsv);
}

start();
