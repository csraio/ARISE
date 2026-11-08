"""Builds ARISE's residual-solvent and water 1H peaks from published tables (dev-time only).

Run with the pinned venv from scripts/install_reference_tools.sh:
    ~/.local/share/arise-reference/bin/python -I scripts/build_solvent_peaks.py

No shift in the output is typed by hand. Three independently typeset tables are read from
their PDF text layers with pdfminer.six:
  gottlieb  Gottlieb, Kotlyar, Nudelman, J. Org. Chem. 1997, 62, 7512 (Table 1)
  fulmer    Fulmer et al., Organometallics 2010, 29, 2176 (Table 1)
  cil       Cambridge Isotope Laboratories, NMR Solvent Data Chart (2010)
A shift is accepted only where at least two sources agree within 0.01 ppm; the accepted
value is the mean of the agreeing values. Fulmer reprints Gottlieb's measurements for the
seven solvents both cover, so agreement there guards the typesetting and the extraction,
not the measurement; CIL measured its own. CIL's water column is "approximate values only"
(its own footnote) and is not used. Each source's column heading is identified with a
PubChem compound by name (ALIASES for headings PubChem does not resolve, or resolves to the
wrong compound), checked against the heading's own formula or -dN deuterium count.

Writes:
  data/solvent_peaks.tsv                  the accepted peaks, compiled into the engine
  tests/data/reference/solvent_sources.tsv every value each source gives, and the verdict

ARISE never imports this script; it only reads the files it writes.
"""

import collections
import datetime
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from decimal import Decimal
from pathlib import Path

from pdfminer.high_level import extract_pages, extract_text
from pdfminer.layout import LAParams, LTChar, LTLine, LTRect, LTTextLine
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "solvent_peaks.tsv"
REPORT = REPO / "tests" / "data" / "reference" / "solvent_sources.tsv"
CACHE = Path.home() / ".cache" / "arise" / "solvents"

# The Pitt server times out over HTTPS from some networks; the pinned sha256 below guards
# the content whichever way it arrives.
SOURCES = {
    "gottlieb": {
        "url": "http://ccc.chem.pitt.edu/wipf/Web/4505.pdf",
        "sha256": "d9e426b396fa417056d3bd60a23f447f6bd77e7616371ccf0d527c073ca0492d",
        "citation": "H. E. Gottlieb, V. Kotlyar, A. Nudelman, J. Org. Chem. 1997, 62, "
                    "7512-7515, Table 1",
    },
    "fulmer": {
        "url": "http://ccc.chem.pitt.edu/wipf/Web/NMR_Impurities.pdf",
        "sha256": "27b0cbfb5a8cfcc53b2e5ca1307e992b8b4c661fb04ea3ded181b8c7063f386c",
        "citation": "G. R. Fulmer et al., Organometallics 2010, 29, 2176-2179, Table 1",
    },
    "cil": {
        "url": "https://chem.umd.edu/sites/default/files/NMR/isotope_solvent.pdf",
        "sha256": "6e8599b123901ad472ca4598782fcbf1da9eaac8a7b3e1cb21c04c7ae2a38e13",
        "citation": "Cambridge Isotope Laboratories, NMR Solvent Data Chart "
                    "(31175_CIL NMR Chart R2, 5/27/10)",
    },
}

# Headings PubChem's name search does not resolve, or resolves to the wrong compound
# ("Trifluoroethanol-d3" gives a hypofluorite), and the name looked up instead. Every
# resolved structure is still checked against the heading's formula or -dN count.
ALIASES = {
    "(CD3)2SO": "dimethyl sulfoxide-d6",
    "C6D5Cl": "chlorobenzene-d5",
    "TFE-d3": "2,2,2-trifluoroethanol-d3",
    "Trifluoroethanol-d3": "2,2,2-trifluoroethanol-d3",
}

TOLERANCE = Decimal("0.01")
NUMBER = re.compile(r"^(\d+\.\d+)[a-z*]?$")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(name):
    source = SOURCES[name]
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{name}.pdf"
    if not path.exists():
        print(f"downloading {source['url']}", file=sys.stderr)
        with urllib.request.urlopen(source["url"], timeout=300) as response:
            path.write_bytes(response.read())
    digest = sha256(path)
    if digest != source["sha256"]:
        sys.exit(f"{name}: sha256 {digest} is not the pinned {source['sha256']}; the file "
                 "changed: check it, then update SOURCES.")
    return path


# ---- words with positions ----------------------------------------------------------------

Word = collections.namedtuple("Word", "text x0 x1 y")


def text_lines(obj):
    if isinstance(obj, LTTextLine):
        yield obj
    elif hasattr(obj, "__iter__"):
        for child in obj:
            yield from text_lines(child)


def words(page):
    """Every word on the page, split at spaces and at gaps wider than 2.5 pt."""
    found = []
    for line in text_lines(page):
        current = []
        for char in line:
            gap = isinstance(char, LTChar) and current and char.x0 - current[-1].x1 > 2.5
            if not isinstance(char, LTChar) or not char.get_text().strip() or gap:
                if current:
                    found.append(current)
                current = [char] if isinstance(char, LTChar) and char.get_text().strip() else []
            else:
                current.append(char)
        if current:
            found.append(current)
    return [Word("".join(c.get_text() for c in w), w[0].x0, w[-1].x1,
                 round((w[0].y0 + w[0].y1) / 2, 1)) for w in found]


def rows(page_words, tolerance=2.0):
    """Words grouped into rows by height on the page, top to bottom, each left to right."""
    ordered = sorted(page_words, key=lambda w: (-w.y, w.x0))
    grouped = []
    for word in ordered:
        if grouped and abs(grouped[-1][0].y - word.y) <= tolerance:
            grouped[-1].append(word)
        else:
            grouped.append([word])
    return [sorted(row, key=lambda w: w.x0) for row in grouped]


def row_text(row):
    return " ".join(w.text for w in row)


def nearest_column(word, columns):
    centre = (word.x0 + word.x1) / 2
    return min(columns, key=lambda c: abs((c.x0 + c.x1) / 2 - centre)).text


def values_by_column(row, columns):
    found = collections.defaultdict(list)
    for word in row:
        match = NUMBER.match(word.text)
        if match:
            found[nearest_column(word, columns)].append(Decimal(match.group(1)))
    return found


def page_text(page_words):
    return " ".join(row_text(r) for r in rows(page_words))


# ---- the three tables --------------------------------------------------------------------

def table_header(page_rows):
    for row in page_rows:
        if [w.text for w in row[:2]] == ["proton", "mult"]:
            return row[2:]
    sys.exit("no 'proton mult' header row")


def table_page(path):
    for page in extract_pages(path, laparams=LAParams(line_margin=0.1, char_margin=1.0)):
        page_words = words(page)
        if any(row_text(r).startswith("proton mult") for r in rows(page_words)):
            return page_words
    sys.exit(f"{path}: no table with a 'proton mult' header (Table 1, 1H NMR data)")


def hdo_footnote(text, pattern):
    match = re.search(pattern, text)
    if not match:
        sys.exit("HDO footnote not found")
    return match


def read_gottlieb(path):
    page_words = table_page(path)
    page_rows = rows(page_words)
    columns = table_header(page_rows)
    residual = water = None
    for row in page_rows:
        text = row_text(row)
        if text.startswith("solvent residual peak"):
            residual = values_by_column(row, columns)
        elif text.startswith("H2O s"):
            water = values_by_column(row, columns)
    note = hdo_footnote(page_text(page_words),
                        r"appears at (\d+\.\d+) and (\d+\.\d+) ppm in acetone and DMSO")
    water["(CD3)2CO"].append(Decimal(note.group(1)))
    water["(CD3)2SO"].append(Decimal(note.group(2)))
    return [w.text for w in columns], residual, water


def read_fulmer(path):
    page_words = table_page(path)
    page_rows = rows(page_words)
    columns = table_header(page_rows)
    residual = collections.defaultdict(list)
    water = None
    in_residual = False
    for row in page_rows:
        text = row_text(row)
        if text.startswith("solvent residual"):
            in_residual = True
        elif text.startswith("water OH s"):
            water = values_by_column(row, columns)
            break
        if in_residual:
            for column, values in values_by_column(row, columns).items():
                residual[column] += values
    note = hdo_footnote(page_text(page_words),
                        r"HDO is also observed in \(CD3\)2SO \((\d+\.\d+) ppm\) and "
                        r"\(CD3\)2CO \((\d+\.\d+) ppm\)")
    water["(CD3)2SO"].append(Decimal(note.group(1)))
    water["(CD3)2CO"].append(Decimal(note.group(2)))
    return [w.text for w in columns], residual, water


def read_cil(path):
    """Rows of the chart lie between full-width rules; the residual 1H column is the first
    '1H Chemical Shift' heading. Values marked (DSS) or (TSP) are referenced to those, not
    TMS, and are left out; so is the approximate HOD column."""
    page = next(iter(extract_pages(path, laparams=LAParams(line_margin=0.1, char_margin=1.0))))
    page_words = words(page)
    rules = sorted({round(o.y0, 1) for o in page if isinstance(o, (LTLine, LTRect))
                    and o.height < 1.5 and o.width > 500}, reverse=True)
    header = next(r for r in rows(page_words) if row_text(r).startswith("1H Chemical Shift"))
    heading = header[0]
    jhd = min((w for w in page_words if w.text == "JHD"), key=lambda w: w.x0)
    rules = [y for y in rules if y < heading.y]
    left = collections.Counter(round(o.x0) for o in page if isinstance(o, (LTLine, LTRect))
                               and o.height < 1.5 and o.width > 500).most_common(1)[0][0]
    rules = [round(o.y0, 1) for o in page if isinstance(o, (LTLine, LTRect))
             and o.height < 1.5 and o.width > 500 and round(o.x0) == left]
    rules = sorted({y for y in rules if y < heading.y}, reverse=True)
    names, residual = [], {}
    for top, bottom in zip(rules, rules[1:]):
        band = [w for w in page_words if bottom < w.y < top]
        label = " ".join(w.text for w in sorted(band, key=lambda w: (-w.y, w.x0))
                         if w.x1 < heading.x0 - 2)
        if not label or "Chemical" in label:
            continue
        name = re.sub(r"-d (\d)", r"-d\1", label)
        cells = sorted((w for w in band if heading.x0 - 2 <= w.x0 < jhd.x0 - 2),
                       key=lambda w: (-w.y, w.x0))
        shifts = []
        for i, word in enumerate(cells):
            match = NUMBER.match(word.text)
            following = cells[i + 1].text if i + 1 < len(cells) else ""
            if match and following not in ("(DSS)", "(TSP)"):
                shifts.append(Decimal(match.group(1)))
        names.append(name)
        residual[name] = shifts
    return names, residual, {}


# ---- identity ----------------------------------------------------------------------------

def pubchem(name):
    """(CID, SMILES) of the one compound PubChem gives for a name, or the reason not."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name, safe="") + "/property/IsomericSMILES/JSON")
    time.sleep(0.25)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            found = json.load(response)["PropertyTable"]["Properties"]
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None, "PubChem does not know the name"
        raise
    if len(found) != 1:
        return None, f"PubChem gives {len(found)} compounds"
    record = found[0]
    return (record["CID"], record.get("IsomericSMILES") or record["SMILES"]), None


def formula_counts(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    counts = collections.Counter()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        counts["D" if symbol == "H" and atom.GetIsotope() == 2 else symbol] += 1
    return counts


def heading_counts(heading):
    """Element counts of a condensed formula heading such as '(CD3)2CO', or None."""
    tokens = re.findall(r"\(|\)\d*|[A-Z][a-z]?\d*", heading)
    if "".join(tokens) != heading:
        return None
    stack = [collections.Counter()]
    for token in tokens:
        if token == "(":
            stack.append(collections.Counter())
        elif token.startswith(")"):
            group = stack.pop()
            times = int(token[1:] or 1)
            for element, n in group.items():
                stack[-1][element] += n * times
        else:
            element = re.match(r"[A-Z][a-z]?", token).group(0)
            stack[-1][element] += int(token[len(element):] or 1)
    counts = stack[0]
    return counts if set(counts) <= {"C", "D", "O", "S", "N", "Cl"} else None


def identify(heading):
    """((CID, SMILES), None) for a column heading, checked against the heading's own
    formula or -dN count (a bare -d is one), or (None, the reason it is not identified)."""
    lookup = ALIASES.get(heading, heading)
    found, reason = pubchem(lookup)
    if found is None:
        return None, reason
    cid, smiles = found
    counts = formula_counts(smiles)
    if counts["H"]:
        return None, f"CID {cid} {smiles} is not fully deuterated"
    expected = heading_counts(heading)
    d_suffix = re.search(r"-d(\d*)$", heading, re.IGNORECASE)
    if expected is not None:
        if expected != +counts:
            return None, f"CID {cid} {smiles} does not have the heading's formula"
    elif d_suffix:
        if int(d_suffix.group(1) or 1) != counts["D"]:
            return None, f"CID {cid} {smiles} has {counts['D']} D, not the heading's"
    else:
        return None, f"neither a formula nor -dN: nothing to check CID {cid} against"
    return (cid, smiles), None


def exchanges_labile_h(smiles):
    """True when the solvent's structure carries deuterium on O, N or S (D2O, CD3OD,
    TFE-d3); ARISE then leaves the suggested exchangeable count open (the chemist's
    decision, 2026-11-01)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return any(atom.GetSymbol() == "H" and atom.GetIsotope() == 2
               and atom.GetNeighbors()[0].GetSymbol() in ("O", "N", "S")
               for atom in mol.GetAtoms())


# ---- agreement ---------------------------------------------------------------------------

def agree(by_source):
    """(accepted shifts, shifts only one source lists), each sorted, from {source: [shifts]}.
    Going up the shift scale, a signal is a group of values from different sources lying
    within TOLERANCE of the group's lowest value (each source's nearest); it is accepted,
    as the group's mean, when at least two sources give it."""
    pending = sorted((v, s) for s, values in by_source.items() for v in values)
    accepted, unmatched = [], []
    while pending:
        low = pending[0][0]
        group = {}
        for value, source in pending:
            if value - low <= TOLERANCE and source not in group:
                group[source] = value
        for source, value in group.items():
            pending.remove((value, source))
        if len(group) >= 2:
            accepted.append(sum(group.values()) / len(group))
        else:
            unmatched += group.values()
    return accepted, unmatched


def printed(values):
    """Shifts as the source prints them."""
    return ",".join(str(v) for v in values) if values else "-"


def fmt(values):
    """Accepted shifts: means to 3 decimals, trailing zeros dropped."""
    return ",".join(f"{v:.3f}".rstrip("0").rstrip(".") for v in values) if values else "-"


# What the sources themselves say about the water peak, quoted in data/solvent_peaks.tsv.
# The build stops if a PDF no longer contains its phrase.
WATER_STATEMENTS = [
    ("gottlieb", "the chemical shift of the residual HDO is very temperature-dependent"),
    ("cil", "may vary with pH"),
]


def water_statements():
    """Header lines quoting WATER_STATEMENTS, each checked against its PDF's text."""
    lines = []
    for name, phrase in WATER_STATEMENTS:
        text = re.sub(r"\s+", " ", extract_text(str(fetch(name))))
        if phrase not in text:
            sys.exit(f"{name}: the PDF no longer says '{phrase}'; update WATER_STATEMENTS")
        lines.append(f"#   {name}: '{phrase}'")
    return lines


def main():
    tables = {
        "gottlieb": read_gottlieb(fetch("gottlieb")),
        "fulmer": read_fulmer(fetch("fulmer")),
        "cil": read_cil(fetch("cil")),
    }
    solvents = {}  # CID -> {"smiles", "headings": {source: heading}}
    unresolved = []
    for source, (headings, _, _) in tables.items():
        for heading in headings:
            found, reason = identify(heading)
            if found is None:
                unresolved.append(f"{source} {heading!r} ({reason})")
                continue
            cid, smiles = found
            entry = solvents.setdefault(cid, {"smiles": smiles, "headings": {}})
            entry["headings"][source] = heading

    report = []
    data = []
    for cid, entry in sorted(solvents.items(), key=lambda kv: kv[0]):
        headings = entry["headings"]
        display = headings.get("fulmer") or headings.get("gottlieb") or headings["cil"]
        accepted = {}
        for kind, index in (("residual", 1), ("water", 2)):
            by_source = {s: tables[s][index].get(h, []) for s, h in headings.items()}
            for source, values in by_source.items():
                report.append("\t".join([display, str(cid), kind, source, headings[source],
                                         printed(sorted(values))]))
            accepted[kind], unmatched = agree(by_source)
            report.append("\t".join([display, str(cid), kind, "accepted", "-",
                                     fmt(accepted[kind])]))
            report.append("\t".join([display, str(cid), kind, "one source only", "-",
                                     printed(sorted(unmatched))]))
        if accepted["residual"]:
            data.append("\t".join([display, str(cid), fmt(accepted["residual"]),
                                   fmt(accepted["water"]),
                                   "yes" if exchanges_labile_h(entry["smiles"]) else "no"]))

    today = datetime.date.today().isoformat()
    sources = [f"# Source {name}: {s['citation']}; {s['url']} sha256={s['sha256']}"
               for name, s in SOURCES.items()]
    common = [f"# Generated by scripts/build_solvent_peaks.py on {today}; do not edit."]
    common += sources
    common += [f"# Aliases looked up on PubChem: " + "; ".join(
        f"{k} -> {v}" for k, v in ALIASES.items())]
    REPORT.write_text("\n".join(common + [
        "# Every residual and water 1H value each source gives, per solvent (PubChem CID),",
        "# what was accepted. Columns: solvent, CID, kind, source (or 'accepted', the means",
        "# of signals two sources give within 0.01 ppm; or 'one source only', values no",
        "# other source confirms), the source's heading, shifts (ppm, as printed; '-' none).",
        "# Water includes each source's HDO footnote.",
    ] + [f"# Not identified: {u}" for u in unresolved] + report) + "\n")
    DATA.write_text("\n".join(common + [
        "# Residual-solvent and water 1H peaks (ppm from TMS) accepted where at least two",
        "# sources agree within 0.01 ppm (mean of the agreeing values). On the water peak the",
        "# sources say (quoted from the PDFs, checked by the builder):",
        *water_statements(),
        "# so a water match means 'possibly water'. A solvent is listed when at least one",
        "# residual peak is accepted; '-' = no",
        "# water value accepted. Signals only one source gives are left out; every value and",
        "# its fate: tests/data/reference/solvent_sources.tsv.",
        "# Columns: solvent, PubChem CID, residual shifts, water shifts, exchanges labile H",
        "# (deuterium on O/N/S in the PubChem structure).",
    ] + data) + "\n")
    print(f"wrote {len(data)} solvents to {DATA}", file=sys.stderr)


if __name__ == "__main__":
    main()
