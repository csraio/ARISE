"""Builds experimental 1H, 13C, DEPT and HSQC reference cases from BMRB (dev-time only).

Run with the pinned venv from scripts/install_reference_tools.sh:
    ~/.local/share/arise-reference/bin/python -I scripts/build_bmrb_reference.py

The cases are chosen by a rule, not by hand: every sample of a BMRB metabolomics entry
(https://bmrb.io, NMR-STAR via the public REST API https://api.bmrb.io/v2) that has
  * a peak list named "1D 1H" and a 1H-13C HSQC peak list (two dimensions, H and C, with
    "HSQC" in the experiment name), both recorded on that sample and both listing peaks
    (some entries give only unassigned spectral transitions);
  * exactly one solute, and exactly one solvent that is one of ARISE's solvents
    (data/solvent_peaks.tsv): the solvent's name, looked up in PubChem, must have the same
    InChI skeleton (formula, connections, hydrogens) as the ARISE solvent's PubChem CID, and
    be deuterated (by name, as "D2O", or by the sample's own "deuterated" labelling, or,
    for a name in arise_solvents.UNLABELLED_DEUTERATED given with no labelling, by
    that rule; the matching is shared in scripts/arise_solvents.py);
  * exactly one compound (chem_comp) whose atoms and bonds, as the entry lists them, give
    RDKit the InChI the entry itself states, and which is in ARISE's scope (net-neutral,
    one fragment, ARISE's elements);
taken in entry-ID order. Every sample left out is counted, by reason, in the file headers.

No shift, assignment or structure is typed by hand: all come from the entry. Atom
assignments are mapped to the heavy atoms of the entry's own structure (an H to the atom
it is bonded to), numbered in the order of the Kekule SMILES written beside them.

Writes tests/data/reference/:
  bmrb_cases.tsv    one row per case: entry, sample, solvent, compound name, Kekule SMILES,
                    and per SMILES atom its BMRB atom name and bonded H count
  bmrb_peaks.tsv    every peak of the case's 1H, 13C, DEPT-135, DEPT-90 and HSQC lists, with
                    the SMILES atoms BMRB assigns it to
  BMRB_NOTICE.md    BMRB's statements on the data, quoted from its pages on the build date,
                    the citation BMRB asks for, and every entry's ID, DOI and authors

ARISE never imports this script; it only reads the files it writes.
"""

import collections
import datetime
import html
import importlib.util
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


def load_solvent_rules():
    """scripts/arise_solvents.py, loaded by path (python -I leaves scripts/ off sys.path)."""
    path = Path(__file__).resolve().parent / "arise_solvents.py"
    spec = importlib.util.spec_from_file_location("arise_solvents", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SOLVENT_RULES = load_solvent_rules()
inchi_skeleton = SOLVENT_RULES.inchi_skeleton

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "tests" / "data" / "reference"
SOLVENTS = REPO / "data" / "solvent_peaks.tsv"
CACHE = Path.home() / ".cache" / "arise" / "bmrb"
PUBCHEM_CACHE = CACHE / "pubchem"
API = "https://api.bmrb.io/v2"
HEADERS = {"Application": "ARISE reference builder"}

ALLOWED = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"}

# BMRB's own statements about its data, quoted in BMRB_NOTICE.md. Each is found on its page
# by the words given here and quoted as the whole sentence the page prints around them; the
# build stops if a page no longer says it.
STATEMENTS = [
    ("https://bmrb.io/bmrb/mission_statement.shtml", "BMRB Mission Statement",
     "worldwide in the public domain"),
    ("https://bmrb.io/metabolomics/data_policy.shtml", "BMRB Metabolomics standards data policy",
     "publicly available free of charge"),
    ("https://bmrb.io/bmrb/aims_and_policies.shtml", "BMRB Aims and Policies",
     "requested to acknowledge BMRB"),
    ("https://bmrb.io/bmrb/aims_and_policies.shtml", "BMRB Aims and Policies",
     "have the author(s)' name(s) attached"),
]
CITATION_URL = "https://bmrb.io/bmrb/citation.shtml"

# Peak lists, by experiment name (BMRB names one experiment several ways).
ONE_D = {"1D 1H": "1H", "1D 13C": "13C", "1D DEPT135": "DEPT135", "1D DEPT90": "DEPT90"}

DROPPED = collections.Counter()
# Samples a rule above admitted, by rule (written into the file headers).
ADMITTED = collections.Counter()
# Optional lists (13C, DEPT-135, DEPT-90) left out of a case because they cannot be read.
LEFT_OUT = collections.Counter()


# ---- fetching ----------------------------------------------------------------------------

def get(url, path=None, delay=0.25):
    """The response body, from the cache when `path` is given and present."""
    if path is not None and path.exists():
        return path.read_bytes()
    time.sleep(delay)
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=300) as response:
        body = response.read()
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
    return body


def entry(entry_id):
    body = get(f"{API}/entry/{entry_id}", CACHE / f"{entry_id}.json")
    return json.loads(body)[entry_id]


def page_text(url):
    """A page's visible text, whitespace collapsed. Never cached: the notice quotes the
    page as it reads on the build date."""
    raw = get(url).decode("utf-8", errors="replace")
    raw = re.sub(r"(?s)<(script|style)\b.*?</\1>", " ", raw)
    raw = re.sub(r"(?s)<!--.*?-->", " ", raw)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", raw))).strip()


def sentence_with(text, words, title):
    """The sentence of `text` containing `words`, as printed (starting after the page's
    title when no full stop comes between them)."""
    at = text.find(words)
    if at < 0:
        return None
    start = text.rfind(". ", 0, at) + 2 if text.rfind(". ", 0, at) >= 0 else 0
    heading = text.rfind(title, 0, at)
    if heading >= 0:
        start = max(start, heading + len(title) + 1)
    end = text.find(". ", at)
    return text[start:(len(text) if end < 0 else end + 1)].strip()


# ---- NMR-STAR access ---------------------------------------------------------------------

def frames(data, category):
    return [sf for sf in data["saveframes"] if sf["category"] == category]


def tags(frame):
    return {tag: value for tag, value in frame["tags"]}


def loop(frame, category):
    for found in frame["loops"]:
        if found["category"] == category:
            return [dict(zip(found["tags"], row)) for row in found["data"]]
    return []


def given(value):
    return value not in (None, "", ".", "?")


# ---- solvents ----------------------------------------------------------------------------

def solvent_of(components, solvents):
    """(ARISE solvent name, None), or (None, the reason the sample's solvent is not one)."""
    listed = [c for c in components if c.get("Type", "").lower() == "solvent"]
    if len(listed) != 1:
        return None, f"{len(listed)} solvents listed"
    name = listed[0].get("Mol_common_name", "")
    label = listed[0].get("Isotopic_labeling")
    solvent, by_rule, reason = SOLVENT_RULES.identify(
        name if given(name) else "", label if given(label) else None, solvents, PUBCHEM_CACHE)
    if by_rule:
        ADMITTED[f'unlabelled "{name}" taken as deuterated (UNLABELLED_DEUTERATED)'] += 1
    return solvent, reason


# ---- structure ---------------------------------------------------------------------------

BOND_ORDERS = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
               "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}


def structure(comp):
    """(RDKit molecule with explicit H, {atom name: index}) from the entry's chem_comp atoms
    and bonds, or (None, reason)."""
    atoms = loop(comp, "_Chem_comp_atom")
    bonds = loop(comp, "_Chem_comp_bond")
    if not atoms or not bonds:
        return None, "no atoms or bonds"
    mol = Chem.RWMol()
    index = {}
    for atom in atoms:
        symbol = atom["Type_symbol"].capitalize()
        if atom["Atom_ID"] in index or Chem.GetPeriodicTable().GetAtomicNumber(symbol) <= 0:
            return None, "bad atom list"
        created = Chem.Atom(symbol)
        charge = atom.get("Charge", "0") if given(atom.get("Charge")) else "0"
        if not re.fullmatch(r"[+-]?\d+", charge):
            return None, "an atom charge that is not a whole number"
        created.SetFormalCharge(int(charge))
        created.SetNoImplicit(True)
        created.SetProp("bmrb_name", atom["Atom_ID"])
        index[atom["Atom_ID"]] = mol.AddAtom(created)
    for bond in bonds:
        order = BOND_ORDERS.get(bond.get("Value_order", "").upper())
        ends = (bond.get("Atom_ID_1"), bond.get("Atom_ID_2"))
        if order is None or not all(end in index for end in ends):
            return None, "bond order not single/double/triple/aromatic, or unknown atom"
        mol.AddBond(index[ends[0]], index[ends[1]], order)
        if order == Chem.BondType.AROMATIC:  # RDKit assigns the Kekule form when sanitizing
            mol.GetBondBetweenAtoms(index[ends[0]], index[ends[1]]).SetIsAromatic(True)
            for end in ends:
                mol.GetAtomWithIdx(index[end]).SetIsAromatic(True)
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # RDKit raises several types; any means the drawing is not valid
        return None, "RDKit rejects the structure"
    return mol.GetMol(), index


def in_scope(mol):
    return ({a.GetSymbol() for a in mol.GetAtoms()} <= ALLOWED
            and sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0
            and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            and len(Chem.GetMolFrags(mol)) == 1
            and any(a.GetSymbol() == "C" for a in mol.GetAtoms()))


def kekule_heavy(mol):
    """(Kekule SMILES of the heavy atoms, the SMILES atoms in order as atoms of `mol`)."""
    heavy = Chem.RemoveHs(mol, sanitize=True)
    Chem.Kekulize(heavy, clearAromaticFlags=True)
    smiles = Chem.MolToSmiles(heavy, kekuleSmiles=True, isomericSmiles=False, canonical=True)
    order = [int(i) for i in re.findall(r"\d+", heavy.GetProp("_smilesAtomOutputOrder"))]
    by_name = {a.GetProp("bmrb_name"): a.GetIdx() for a in mol.GetAtoms()}
    return smiles, [mol.GetAtomWithIdx(by_name[heavy.GetAtomWithIdx(i).GetProp("bmrb_name")])
                    for i in order]


# ---- peak lists --------------------------------------------------------------------------

def peak_list_kind(frame):
    info = tags(frame)
    name = info.get("Experiment_name", "")
    dims = sorted(d.get("Atom_type", "") for d in loop(frame, "_Spectral_dim"))
    if name in ONE_D and len(dims) == 1:
        return ONE_D[name]
    if "HSQC" in name and dims == ["C", "H"]:
        return "HSQC"
    return None


def peaks(frame, heavy_of, kind):
    """Rows (h shift, c shift, phase, integral, h atoms, c atoms) of one list, or None when
    an assignment names an atom the structure lacks."""
    dims = {d["ID"]: d["Atom_type"] for d in loop(frame, "_Spectral_dim")}

    def dim(row):
        """The row's dimension (atom type); a 1D list may omit the column."""
        listed = row.get("Spectral_dim_ID")
        if given(listed):
            return dims.get(listed)
        return next(iter(dims.values())) if len(dims) == 1 else None

    shift = collections.defaultdict(dict)
    phase = {}
    for row in loop(frame, "_Peak_char"):
        if dim(row) is None or not re.fullmatch(r"-?\d+(\.\d+)?",
                                                row.get("Chem_shift_val", "")):
            return None
        shift[row["Peak_ID"]][dim(row)] = row["Chem_shift_val"]
        if given(row.get("Phase_val", ".")):
            phase[row["Peak_ID"]] = {"positive": "positive", "negative": "negative"}.get(
                row["Phase_val"].lower(), "?")
    integral = {row["Peak_ID"]: row["Intensity_val"]
                for row in loop(frame, "_Peak_general_char")
                if row.get("Measurement_method", "").lower() == "integration"
                and given(row.get("Intensity_val"))}
    atoms = collections.defaultdict(lambda: {"H": set(), "C": set()})
    for row in loop(frame, "_Assigned_peak_chem_shift"):
        name = row.get("Atom_ID")
        if not given(name):
            continue
        if name not in heavy_of:
            return None
        if dim(row) is None:
            return None
        atoms[row["Peak_ID"]][dim(row)].add(heavy_of[name])
    rows = []
    for peak_id, by_dim in shift.items():
        if kind == "DEPT135" and phase.get(peak_id) not in ("positive", "negative"):
            return None
        if kind == "HSQC" and not ("H" in by_dim and "C" in by_dim):
            return None  # a cross peak needs both shifts
        rows.append((by_dim.get("H", "-"), by_dim.get("C", "-"), phase.get(peak_id, "-"),
                     integral.get(peak_id, "-"), atoms[peak_id]["H"], atoms[peak_id]["C"]))
    return sorted(rows, key=lambda r: [float(v) if v != "-" else 0.0 for v in r[:2]])


def numbers(values):
    return ",".join(str(v) for v in sorted(values)) or "-"


# ---- cases -------------------------------------------------------------------------------

def cases_of(entry_id, data, solvents):
    """The entry's cases, as (case row, peak rows), counting each sample left out."""
    lists = collections.defaultdict(dict)
    for frame in frames(data, "spectral_peak_list"):
        kind = peak_list_kind(frame)
        if kind is None:
            continue
        sample = tags(frame).get("Sample_ID")
        if kind in lists[sample]:
            lists[sample][kind] = None  # two lists of one kind on one sample: ambiguous
        else:
            lists[sample][kind] = frame
    samples = [s for s, k in sorted(lists.items()) if k.get("1H") and k.get("HSQC")]
    if not samples:
        DROPPED["no sample with one 1D 1H and one HSQC list"] += 1
        return []
    comps = frames(data, "chem_comp")
    if len(comps) != 1:
        DROPPED["not exactly one compound"] += len(samples)
        return []
    info = tags(comps[0])
    mol, index = structure(comps[0])
    if mol is None:
        DROPPED[f"structure: {index}"] += len(samples)
        return []
    listed = info.get("InChI_code", "")
    if not given(listed) or inchi_skeleton(Chem.MolToInchi(mol)) != inchi_skeleton(listed):
        DROPPED["structure does not give the entry's InChI"] += len(samples)
        return []
    if not in_scope(mol):
        DROPPED["compound outside ARISE's scope"] += len(samples)
        return []
    smiles, order = kekule_heavy(mol)
    position = {atom.GetIdx(): n for n, atom in enumerate(order)}
    heavy_of = {}
    for atom in mol.GetAtoms():
        name = atom.GetProp("bmrb_name")
        if atom.GetSymbol() == "H":
            if atom.GetDegree() != 1:
                DROPPED["an H bonded to more than one atom"] += len(samples)
                return []
            heavy_of[name] = position[atom.GetNeighbors()[0].GetIdx()]
        else:
            heavy_of[name] = position[atom.GetIdx()]
    names = " ".join(atom.GetProp("bmrb_name") for atom in order)
    hydrogens = " ".join(str(sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == "H"))
                         for atom in order)
    components = collections.defaultdict(list)
    for frame in frames(data, "sample"):
        for row in loop(frame, "_Sample_component"):
            components[row.get("Sample_ID", tags(frame).get("ID"))].append(row)
    found = []
    for sample in samples:
        solutes = [c for c in components[sample] if c.get("Type", "").lower() == "solute"]
        if len(solutes) != 1:
            DROPPED["not exactly one solute"] += 1
            continue
        solvent, reason = solvent_of(components[sample], solvents)
        if solvent is None:
            DROPPED[reason] += 1
            continue
        peak_rows = []
        for kind in ("1H", "13C", "DEPT135", "DEPT90", "HSQC"):
            frame = lists[sample].get(kind)
            if frame is None:
                continue
            rows = peaks(frame, heavy_of, kind)
            if rows == [] and kind in ("1H", "HSQC"):
                DROPPED[f"{kind} list has no peaks (only spectral transitions)"] += 1
                peak_rows = None
                break
            if rows == []:
                LEFT_OUT[f"{kind} (no peaks)"] += 1
                continue
            if rows is None and kind in ("1H", "HSQC"):
                DROPPED[f"{kind} list unusable (unknown atom, dimension, shift or phase)"] += 1
                peak_rows = None
                break
            if rows is None:
                LEFT_OUT[kind] += 1
                continue
            peak_rows += [(kind, n + 1, *row) for n, row in enumerate(rows)]
        if peak_rows is None:
            continue
        found.append(((entry_id, sample, solvent, info.get("Name", "-"), smiles, names,
                       hydrogens), peak_rows))
    return found


# ---- notice ------------------------------------------------------------------------------

def notice(today, authors):
    lines = [
        "# BMRB data in ARISE's tests",
        "",
        "`bmrb_cases.tsv` and `bmrb_peaks.tsv` contain data from the Biological Magnetic "
        "Resonance Data Bank (BMRB, https://bmrb.io), metabolomics entries, read through "
        "BMRB's public API by `scripts/build_bmrb_reference.py`. The data are used only to "
        "test ARISE and are not part of the ARISE engine.",
        "",
        f"BMRB publishes no formal licence. Its own statements, quoted from its pages on "
        f"{today} (generated by the builder, which stops if a page no longer says this):",
        "",
    ]
    for url, title, words in STATEMENTS:
        quote = sentence_with(page_text(url), words, title)
        if quote is None:
            sys.exit(f"{url} no longer contains '{words}': check BMRB's terms, then update "
                     "STATEMENTS")
        lines += [f"> {quote}", ">", f"> — {title}, {url}", ""]
    text = page_text(CITATION_URL)
    start, end = text.find("Citation "), text.find(" Other BMRB Publications")
    if start < 0 or end < start:
        sys.exit(f"{CITATION_URL}: the citation section was not found")
    lines += [
        f"The citation BMRB asks for ({CITATION_URL}):",
        "",
        f"> {text[start + len('Citation '):end].strip()}",
        "",
        "## Entries",
        "",
        "Each entry's own authors and DOI, as BMRB lists them.",
        "",
        "| Entry | DOI | Authors |",
        "|---|---|---|",
    ]
    lines += [f"| {e} | {doi} | {names} |" for e, (doi, names) in sorted(authors.items())]
    return "\n".join(lines) + "\n"


def entry_authors(data):
    info = tags(frames(data, "entry_information")[0])
    names = []
    for row in loop(frames(data, "entry_information")[0], "_Entry_author"):
        parts = [row.get("Given_name"), row.get("Middle_initials"), row.get("Family_name")]
        names.append(" ".join(p for p in parts if given(p)))
    return info.get("DOI", "-"), "; ".join(names) or "-"


# ---- main --------------------------------------------------------------------------------

def tsv(value):
    text = str(value)
    if "\t" in text or "\n" in text:
        sys.exit(f"tab or newline in a value: {text!r}")
    return text


def main():
    today = datetime.date.today().isoformat()
    solvents = SOLVENT_RULES.arise_solvents(SOLVENTS, PUBCHEM_CACHE)
    names = json.loads(get(f"{API}/search/get_all_values_for_tag/"
                           "Spectral_peak_list.Experiment_name?database=metabolomics"))
    candidates = sorted(e for e, listed in names.items()
                        if "1D 1H" in listed and any("HSQC" in n for n in listed))
    cases, peak_rows, authors = [], [], {}
    for entry_id in candidates:
        data = entry(entry_id)
        found = cases_of(entry_id, data, solvents)
        for case, rows in found:
            cases.append(case)
            peak_rows += [(case[0], case[1], *row) for row in rows]
        if found:
            authors[entry_id] = entry_authors(data)
    dropped = [f"#   {n} {reason}" for reason, n in sorted(DROPPED.items())]
    dropped.append("# Optional lists left out of a case (unknown atom, dimension, shift or "
                   "phase): " + (", ".join(f"{kind} {n}" for kind, n in sorted(LEFT_OUT.items()))
                                 or "none"))
    head = [
        f"# Generated by scripts/build_bmrb_reference.py on {today}; do not edit.",
        f"# RDKit {rdBase.rdkitVersion}. Source: BMRB metabolomics entries via {API} "
        "(see BMRB_NOTICE.md).",
        f"# {len(candidates)} entries list a 1D 1H and an HSQC; {len(cases)} samples of "
        f"{len(authors)} entries are cases. Samples left out, by reason:",
        *dropped,
    ]
    head += [f"# Admitted by rule: {n} {rule}" for rule, n in sorted(ADMITTED.items())]
    by_solvent = collections.Counter(case[2] for case in cases)
    head.append("# Cases by solvent: " + ", ".join(f"{s} {n}"
                                                   for s, n in sorted(by_solvent.items())))
    case_lines = head + [
        "# Columns: entry, sample, solvent (ARISE name), compound, Kekule SMILES (heavy "
        "atoms), BMRB atom name per SMILES atom, bonded H count per SMILES atom.",
        *("\t".join(tsv(v) for v in case) for case in cases),
    ]
    (OUT / "bmrb_cases.tsv").write_text("\n".join(case_lines) + "\n")
    peak_lines = head + [
        "# Columns: entry, sample, list (1H|13C|DEPT135|DEPT90|HSQC), peak (1-based, by "
        "shift), 1H shift, 13C shift, phase (positive|negative), integral, SMILES atoms assigned in the 1H "
        "dimension (an H's bonded atom), SMILES atoms assigned in the 13C dimension; '-' "
        "when absent.",
        *("\t".join(tsv(v) for v in (*row[:8], numbers(row[8]), numbers(row[9])))
          for row in peak_rows),
    ]
    (OUT / "bmrb_peaks.tsv").write_text("\n".join(peak_lines) + "\n")
    (OUT / "BMRB_NOTICE.md").write_text(notice(today, authors))
    for line in head:
        print(line, file=sys.stderr)


if __name__ == "__main__":
    main()
