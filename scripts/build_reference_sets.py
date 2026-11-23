"""Builds independent reference data for ARISE's chem tests (dev-time only).

Run with the pinned venv from scripts/install_reference_tools.sh:
    ~/.local/share/arise-reference/bin/python scripts/build_reference_sets.py

Writes tests/data/reference/:
  nmrshiftdb2_13c_sample.tsv   Experimentally assigned 13C peaks (nmrshiftdb2) for a
                               deterministic sample of small molecules: Kekule SMILES and,
                               per atom, the assigned peak ('.' non-carbon, '-' carbon not
                               listed; databases list one atom per environment).
  NMRSHIFTDB2_LICENSE.txt      The nmrshiftdb2 Database License, required with the extract.
  pubchem_rdkit_oracle.tsv     For every cleaned PubChem candidate of C4H8O2, C8H8O2 and
                               C10H12O2: per-atom H counts, symmetry classes and the counts
                               of RDKit's published functional-group fragments (fr_*),
                               computed by RDKit, an implementation independent of ARISE's.
  chebi_classes.tsv            A deterministic sample of ChEBI 3-star molecules with the
                               expert-curated classes (is_a ancestors) they belong to, among
                               the 27 classes of docs/chem_group_definitions.md.
  pubchem_<formula>.csv        The raw PubChem responses the oracle was built from.

No functional-group definition in this script is written by ARISE's authors: the groups
come from RDKit's fragment library and ChEBI's ontology.

ARISE never imports this script or RDKit; tests only read the files it writes.
"""

import argparse
import collections
import datetime
import gzip
import hashlib
import importlib.util
import json
import re
import statistics
import shutil
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import Fragments, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "tests" / "data" / "reference"
NMRSHIFTDB_URL = (
    "https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2.nmredata.sd/download"
)
LICENSE_URL = "https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/nmrshiftdb2datalicense.txt"
PUBCHEM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{}/property/"
    "ConnectivitySMILES/CSV"
)
CHEBI_SDF_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite_3_stars.sdf.gz"
CHEBI_OBO_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo.gz"
CHEBI_STRUCTURES_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/structures.tsv.gz"
PUBCHEM_INCHIKEY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/InChI,InChIKey/JSON"
)
ALLOWED = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"}
SAMPLE_SIZE = 400
MAX_HEAVY_ATOMS = 25

# RDKit fragment functions (rdkit.Chem.Fragments, defined in RDKit's
# FragmentDescriptors.csv) recorded for the oracle. Which ARISE groups each one is compared
# with is decided in tests/chem/reference_test.cpp, where the definitions agree.
FRAGMENTS = ["fr_ketone", "fr_aldehyde", "fr_COO2", "fr_ester", "fr_Al_OH", "fr_phenol",
             "fr_nitrile", "fr_urea", "fr_benzene", "fr_pyridine", "fr_furan", "fr_thiophene",
             "fr_methoxy"]

DEFINITIONS = REPO / "docs" / "chem_group_definitions.md"


def definition_classes():
    """The ChEBI classes of the signed-off groups table, in row order (a row may name several,
    "CHEBI:<id> <name>; CHEBI:<id> <name>", e.g. methoxy). Read from the doc so the IDs are
    typed in one place only (tests compare the C++ table with it)."""
    classes = []
    rows = 0
    for line in DEFINITIONS.read_text(encoding="utf-8").splitlines():
        cells = line.split("|")
        if len(cells) != 9 or not re.fullmatch(r"\s*\d+\s*", cells[1]):
            continue  # not a data row of the 7-column groups table
        ids = re.findall(r"CHEBI:(\d+) ", cells[6])
        if not ids:
            continue
        rows += 1
        if int(cells[1]) != rows:
            sys.exit(f"{DEFINITIONS}: groups table row {cells[1].strip()} out of order")
        classes.extend(int(i) for i in ids)
    if not classes:
        sys.exit(f"{DEFINITIONS}: no groups table found")
    return classes


CHEBI_CLASSES = definition_classes()
# ChEBI's ring classes (furans, pyridines, ...) are skeleton classes that include
# non-aromatic rings (butenolides, dihydropyridines). ARISE's ring groups mean aromatic
# rings, so each ChEBI row also records RDKit's aromatic ring counts, and a ChEBI ring
# class is checked only where RDKit sees the ring (user decision 2026-10-31). Pyrroles have
# no RDKit fragment; the pattern named "pyrrole_ring" is written like RDKit's fr_furan
# pattern, and known_differences.tsv ties the name to the PubChem compound it encodes
# (verify_ring_patterns).
CHEBI_RINGS = [("fr_benzene", None), ("fr_pyridine", None), ("fr_furan", None),
               ("fr_thiophene", None), ("pyrrole_ring", Chem.MolFromSmarts("n1cccc1"))]
CHEBI_PER_CLASS = 100
CHEBI_CHARGED_PER_CLASS = 10
CHEBI_MAX_HEAVY_ATOMS = 40


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, path):
    print(f"downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=900) as response, open(path, "wb") as out:
        while block := response.read(1 << 20):
            out.write(block)


def refresh_license(path):
    """Re-downloads the nmrshiftdb2 licence text to `path`. When the site does not answer,
    the copy already at `path` is kept with a warning, so a rebuild does not depend on the
    site being up; with no copy at all, the build stops."""
    partial = path.with_name(path.name + ".partial")
    try:
        download(LICENSE_URL, partial)
    except OSError as error:  # URLError, timeouts and HTTP errors are all OSErrors
        if partial.exists():
            partial.unlink()
        if not path.exists():
            sys.exit(f"cannot download the nmrshiftdb2 licence ({error}) and no copy exists "
                     f"at {path}")
        print(f"WARNING: cannot download the nmrshiftdb2 licence ({error}); keeping the copy "
              f"already at {path}", file=sys.stderr)
        return
    partial.replace(path)


def header(lines, sources):
    today = datetime.date.today().isoformat()
    head = [
        "# Generated by scripts/build_reference_sets.py on " + today + "; do not edit.",
        f"# RDKit {rdBase.rdkitVersion}.",
    ]
    head += [f"# Source: {source}" for source in sources]
    return "\n".join(head + lines) + "\n"


# ---- record integrity -------------------------------------------------------------------

# Records dropped because the structure RDKit reads from the record does not give the
# InChI the database itself lists for it, per source (written into each file's header).
DROPPED = collections.Counter()


def inchi_skeleton(inchi):
    """The formula, connection and hydrogen layers of an InChI: what ARISE uses (it ignores
    stereo, and nmrshiftdb2 adds an isotope layer of all-natural '+0' entries)."""
    layers = inchi.strip().split("/")
    kept = layers[:2]
    for layer in layers[2:]:
        if layer[:1] not in ("c", "h"):
            break  # main layer ends; later /h layers belong to the isotopic layer
        kept.append(layer)
    return "/".join(kept)


def same_structure(mol, listed_inchi):
    """Whether RDKit's InChI of `mol` agrees with the database's own InChI."""
    return bool(listed_inchi) and inchi_skeleton(Chem.MolToInchi(mol)) == inchi_skeleton(
        listed_inchi)


# ---- nmrshiftdb2 ------------------------------------------------------------------------


def parse_13c_labels(block):
    return set(re.findall(r"L=([^\s,\\]+)", block))


def parse_assignment(block):
    """label -> (shift, [molfile atom numbers])"""
    result = {}
    for raw in block.splitlines():
        line = raw.strip().rstrip("\\").strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 3:
            continue
        try:
            result[parts[0]] = (float(parts[1]), [int(p) for p in parts[2:]])
        except ValueError:
            return None
    return result


def kekule_smiles_with_order(mol):
    """Heavy-atom Kekule SMILES and, per SMILES atom, the molfile atom number."""
    for atom in mol.GetAtoms():
        atom.SetIntProp("molfile_number", atom.GetIdx() + 1)
    heavy = Chem.RemoveHs(mol)
    Chem.Kekulize(heavy, clearAromaticFlags=True)
    smiles = Chem.MolToSmiles(heavy, kekuleSmiles=True, isomericSmiles=False, canonical=True)
    order = [int(i) for i in re.findall(r"\d+", heavy.GetProp("_smilesAtomOutputOrder"))]
    return smiles, [heavy.GetAtomWithIdx(i) for i in order]


def nmrshiftdb_record(mol):
    if not mol.HasProp("NMREDATA_1D_13C") or not mol.HasProp("NMREDATA_ASSIGNMENT"):
        return None
    listed = mol.GetProp("NMREDATA_INCHI") if mol.HasProp("NMREDATA_INCHI") else ""
    if not same_structure(mol, listed):
        DROPPED["nmrshiftdb2"] += 1
        return None
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    heavy_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    if not elements <= ALLOWED or "C" not in elements or heavy_count > MAX_HEAVY_ATOMS:
        return None
    if net_charge(mol) or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    labels = parse_13c_labels(mol.GetProp("NMREDATA_1D_13C"))
    assignment = parse_assignment(mol.GetProp("NMREDATA_ASSIGNMENT"))
    if not labels or assignment is None:
        return None
    peaks = sorted((assignment[l][0], l) for l in labels if l in assignment)
    if len(peaks) < 2:
        return None
    peak_of_atom = {}
    for index, (_, label) in enumerate(peaks):
        for number in assignment[label][1]:
            if number in peak_of_atom or number > mol.GetNumAtoms():
                return None  # one atom on two peaks, or a bad index: unusable record
            if mol.GetAtomWithIdx(number - 1).GetSymbol() != "C":
                return None
            peak_of_atom[number] = index
    smiles, atoms = kekule_smiles_with_order(mol)
    marks = []
    shifts = []
    for atom in atoms:
        if atom.GetSymbol() != "C":
            marks.append(".")
            shifts.append(".")
            continue
        peak = peak_of_atom.get(atom.GetIntProp("molfile_number"))
        marks.append("-" if peak is None else str(peak))
        shifts.append("-" if peak is None else f"{peaks[peak][0]:g}")
    return smiles, " ".join(marks), " ".join(shifts)


def build_nmrshiftdb(sd_path):
    records = []
    supplier = Chem.ForwardSDMolSupplier(str(sd_path), removeHs=False, sanitize=True)
    for mol in supplier:
        if mol is None or not mol.HasProp("NMREDATA_ID"):
            continue
        match = re.search(r"DB_ID=(\d+)", mol.GetProp("NMREDATA_ID"))
        if not match:
            continue
        try:
            record = nmrshiftdb_record(mol)
        except Exception:  # an RDKit failure on one record must not stop the build
            record = None
        if record is not None:
            records.append((int(match.group(1)), *record))
    records.sort()
    step = max(1, len(records) // SAMPLE_SIZE)
    sample = records[::step][:SAMPLE_SIZE]
    print(f"nmrshiftdb2: {len(records)} usable 13C records, sampled {len(sample)}",
          file=sys.stderr)
    return ([f"{db_id}\t{smiles}\t{marks}" for db_id, smiles, marks, _ in sample],
            [f"{db_id}\t{smiles}\t{shifts}" for db_id, smiles, _, shifts in records])


# ---- exchangeable protons (secondary cases for hsqc) ------------------------------------

def load_solvent_rules():
    """scripts/arise_solvents.py, loaded by path (python -I leaves scripts/ off sys.path)."""
    path = Path(__file__).resolve().parent / "arise_solvents.py"
    spec = importlib.util.spec_from_file_location("arise_solvents", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SOLVENT_RULES = load_solvent_rules()
SOLVENT_TABLE = REPO / "data" / "solvent_peaks.tsv"
PUBCHEM_CACHE = Path.home() / ".cache" / "arise" / "pubchem"
EXCHANGING = {"O", "N", "S", "Se"}  # H on these is exchangeable (chem::exchangeable_hydrogens)
EXCHANGEABLE_DROPPED = collections.Counter()
EXCHANGEABLE_ADMITTED = collections.Counter()


def record_solvent(text, solvents):
    """The ARISE solvent an nmrshiftdb2 NMREDATA_SOLVENT text names, or None (counted).
    The whole text is looked up first, then each comma-separated part ("DMSO-D6, C2D6SO"
    gives a name and a formula). Parts PubChem does not know are skipped; at least one part
    must name an ARISE solvent, every part that names a compound must name that same one,
    and a part naming anything else (another or an undeuterated compound) leaves the record
    out."""
    text = text.strip().splitlines()[0].rstrip("\\").strip() if text.strip() else ""
    if not text or text.lower() in ("unreported", "unknown"):
        EXCHANGEABLE_DROPPED["no solvent reported"] += 1
        return None
    unknown = "solvent name not one PubChem compound"
    for names in ([text], [part.strip() for part in text.split(",") if part.strip()]):
        found = set()
        rejected = False
        for name in names:
            solvent, by_rule, reason = SOLVENT_RULES.identify(name, None, solvents, PUBCHEM_CACHE)
            if solvent is not None:
                found.add((solvent, by_rule, name))
            elif reason != unknown:
                rejected = True
        if rejected:
            break
        if len({solvent for solvent, _, _ in found}) == 1:
            solvent, by_rule, name = sorted(found)[0]
            if by_rule:
                EXCHANGEABLE_ADMITTED[f'unlabelled "{name}" taken as deuterated'] += 1
            return solvent
        if len(found) > 1:
            EXCHANGEABLE_DROPPED["solvent mixture"] += 1
            return None
    EXCHANGEABLE_DROPPED["solvent not one of ARISE's (or not deuterated)"] += 1
    return None


def proton_peaks(mol):
    """[(shift, kind, [H-bearing atom indices])] of the record's 1H spectrum, kind "C" for H
    on carbon and "X" for H on O/N/S/Se; None (counted) when a peak is unusable."""
    assignment = parse_assignment(mol.GetProp("NMREDATA_ASSIGNMENT"))
    if assignment is None:
        EXCHANGEABLE_DROPPED["unreadable assignment"] += 1
        return None
    peaks = []
    seen = set()
    for raw in mol.GetProp("NMREDATA_1D_1H").splitlines():
        match = re.match(r"\s*(-?\d+(?:\.\d+)?)\s*,.*\bL=([^\s,\\]+)", raw)
        if not match or match.group(2) in seen:
            continue  # some records repeat a peak's line under the same label: one peak
        label = match.group(2)
        seen.add(label)
        if label not in assignment:
            EXCHANGEABLE_DROPPED["1H peak without an assignment"] += 1
            return None
        bearers = set()
        for number in assignment[label][1]:
            if not 0 < number <= mol.GetNumAtoms():
                EXCHANGEABLE_DROPPED["assignment names a missing atom"] += 1
                return None
            atom = mol.GetAtomWithIdx(number - 1)
            if atom.GetSymbol() == "H":
                if atom.GetDegree() != 1:
                    EXCHANGEABLE_DROPPED["H bonded to more than one atom"] += 1
                    return None
                atom = atom.GetNeighbors()[0]
            bearers.add(atom.GetIdx())
        elements = {mol.GetAtomWithIdx(i).GetSymbol() for i in bearers}
        if elements == {"C"}:
            kind = "C"
        elif elements and elements <= EXCHANGING:
            kind = "X"
        else:
            EXCHANGEABLE_DROPPED["1H peak on H of several kinds, or on P/Si/B"] += 1
            return None
        peaks.append((float(match.group(1)), kind, sorted(bearers)))
    return peaks


def exchangeable_case(mol, solvents):
    """(solvent, SMILES, peaks text, HSQC text) of one record, or None (counted)."""
    for prop in ("NMREDATA_1D_1H", "NMREDATA_1D_13C", "NMREDATA_ASSIGNMENT", "NMREDATA_INCHI"):
        if not mol.HasProp(prop):
            return None
    if not in_scope(mol) or not same_structure(mol, mol.GetProp("NMREDATA_INCHI")):
        return None
    peaks = proton_peaks(mol)
    if not peaks or not any(kind == "X" for _, kind, _ in peaks):
        return None
    solvent = record_solvent(mol.GetProp("NMREDATA_SOLVENT")
                             if mol.HasProp("NMREDATA_SOLVENT") else "", solvents)
    if solvent is None:
        return None
    carbon = carbon_shifts(mol)
    hsqc = set()
    for shift, kind, bearers in peaks:
        if kind != "C":
            continue
        for atom in bearers:
            if atom not in carbon:
                EXCHANGEABLE_DROPPED["C-H peak whose carbon has no assigned 13C shift"] += 1
                return None
            hsqc.add((shift, carbon[atom]))
    if not hsqc:
        EXCHANGEABLE_DROPPED["no H on carbon (its HSQC would be empty)"] += 1
        return None
    smiles, _ = kekule_smiles_with_order(mol)
    peak_text = " ".join(f"{shift:g}:{kind}" for shift, kind, _ in sorted(peaks))
    hsqc_text = " ".join(f"{h:g}:{c:g}" for h, c in sorted(hsqc))
    return solvent, smiles, peak_text, hsqc_text


def build_exchangeable_cases(sd_path):
    """Lines of nmrshiftdb2_exchangeable.tsv: every InChI-checked, in-scope record with an
    assigned 1H peak on H bonded to O, N, S or Se, an assigned 13C spectrum and a solvent
    that is one of ARISE's. Secondary evidence: the HSQC is derived from the assignments."""
    solvents = SOLVENT_RULES.arise_solvents(SOLVENT_TABLE, PUBCHEM_CACHE)
    rows = []
    for mol in Chem.ForwardSDMolSupplier(str(sd_path), removeHs=False, sanitize=True):
        if mol is None or not mol.HasProp("NMREDATA_ID"):
            continue
        match = re.search(r"DB_ID=(\d+)", mol.GetProp("NMREDATA_ID"))
        if not match:
            continue
        try:
            case = exchangeable_case(mol, solvents)
        except Exception:  # an RDKit failure on one record must not stop the build
            EXCHANGEABLE_DROPPED["RDKit failure"] += 1
            case = None
        if case is not None:
            rows.append((int(match.group(1)), *case))
    rows.sort()
    by_solvent = collections.Counter(row[1] for row in rows)
    print(f"nmrshiftdb2: {len(rows)} exchangeable-proton cases", file=sys.stderr)
    head = [f"# Records with an OH/NH/SH/SeH 1H peak left out, by reason (records with no such",
            "# peak, outside the scope or failing the InChI check are not counted):"]
    head += [f"#   {n} {reason}" for reason, n in sorted(EXCHANGEABLE_DROPPED.items())]
    head += [f"# Admitted by rule: {n} {rule} (arise_solvents.UNLABELLED_DEUTERATED)"
             for rule, n in sorted(EXCHANGEABLE_ADMITTED.items())]
    head.append("# Cases by solvent: " + ", ".join(f"{s} {n}" for s, n in sorted(by_solvent.items())))
    return head + ["\t".join(str(v) for v in row) for row in rows]


# Measured, not remembered: the C=N carbon of oximes vs oxime ethers (the user's rule for
# a new group: only if its shift differs appreciably; docs/chem_group_definitions.md), with
# the classes stated as SMARTS here so the measurement is reproducible.
CN_CLASSES = [("oxime", Chem.MolFromSmarts("[CX3]=[NX2][OX2H1]")),
              ("oxime-ether", Chem.MolFromSmarts("[CX3]=[NX2][OX2][#6]"))]


def spectrum_shifts(block):
    """Every peak shift of an NMREDATA_1D_* block (lines 'shift, L=label')."""
    shifts = []
    for raw in block.splitlines():
        match = re.match(r"\s*(-?\d+(?:\.\d+)?)\s*,", raw)
        if match:
            shifts.append(float(match.group(1)))
    return shifts


def carbon_shifts(mol):
    """molfile atom index (0-based) -> assigned 13C shift."""
    assignment = parse_assignment(mol.GetProp("NMREDATA_ASSIGNMENT")) or {}
    labels = parse_13c_labels(mol.GetProp("NMREDATA_1D_13C"))
    shifts = {}
    for label in labels & set(assignment):
        shift, atoms = assignment[label]
        for number in atoms:
            if 0 < number <= mol.GetNumAtoms() and mol.GetAtomWithIdx(number - 1).GetSymbol() == "C":
                shifts[number - 1] = shift
    return shifts


def in_scope(mol):
    """Net-neutral, single-fragment molecules of ARISE's elements (the scope of every
    sample); charge-separated forms such as nitro groups and N-oxides are included."""
    return ({a.GetSymbol() for a in mol.GetAtoms()} <= ALLOWED and net_charge(mol) == 0
            and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            and len(Chem.GetMolFrags(mol)) == 1)


def measure_nmrshiftdb(sd_path):
    """Measurement lines for measurements.tsv: the range of every listed 1H and 13C shift
    of in-scope molecules (the plausibility windows of io/peak_csv.hpp and io/carbon_csv.hpp
    must contain it), and the C=N carbon shifts of CN_CLASSES. Records failing the InChI
    check are skipped."""
    extremes = {}
    cn = collections.defaultdict(list)
    for mol in Chem.ForwardSDMolSupplier(str(sd_path), removeHs=False, sanitize=True):
        if mol is None or not mol.HasProp("NMREDATA_INCHI"):
            continue
        try:
            if not in_scope(mol) or not same_structure(mol, mol.GetProp("NMREDATA_INCHI")):
                continue
        except Exception:  # an RDKit failure on one record must not stop the build
            continue
        for nucleus in ("1H", "13C"):
            prop = f"NMREDATA_1D_{nucleus}"
            if mol.HasProp(prop):
                for shift in spectrum_shifts(mol.GetProp(prop)):
                    low, high = extremes.get(nucleus, (shift, shift))
                    extremes[nucleus] = (min(low, shift), max(high, shift))
        if mol.HasProp("NMREDATA_1D_13C") and mol.HasProp("NMREDATA_ASSIGNMENT"):
            shifts = carbon_shifts(mol)
            for atom in mol.GetAtoms():
                atom.SetIntProp("molfile_index", atom.GetIdx())
            heavy = Chem.RemoveHs(mol)  # SMARTS H counts need implicit hydrogens
            for name, pattern in CN_CLASSES:
                for match in heavy.GetSubstructMatches(pattern):
                    carbon = heavy.GetAtomWithIdx(match[0]).GetIntProp("molfile_index")
                    if carbon in shifts:
                        cn[name].append(shifts[carbon])
    lines = [f"shift-range\t{nucleus}\t{low:g}\t{high:g}"
             for nucleus, (low, high) in sorted(extremes.items())]
    for name, _ in CN_CLASSES:
        values = sorted(cn[name])
        cuts = statistics.quantiles(values, n=20, method="inclusive")
        lines.append(f"cn-shift\t{name}\t{len(values)}\t{statistics.median(values):.1f}"
                     f"\t{cuts[0]:.1f}\t{cuts[-1]:.1f}")
    return lines


# Cut points of the hybridization shift measurement, in thousandths (statistics.quantiles
# n=1000, inclusive): 0.1st, 1st, 99th and 99.9th percentile. Windows use 0.1-99.9 (the
# chemist's decision, 2026-11-22).
HYBRIDIZATION_CUTS = (1, 10, 990, 999)


def measure_hybridization_shifts(extract_lines):
    """Measurement lines for measurements.tsv: the 13C shifts of sp, sp2 and sp3 carbons
    (RDKit's hybridization; uncharged carbons) over the extract lines with even nmrshiftdb2
    IDs, one value per carbon per line. An independent check of tools/hybridization_shifts,
    which counts bonds with ARISE's own reader."""
    values = collections.defaultdict(list)
    molecules = collections.defaultdict(set)
    names = {1: "sp", 2: "sp2", 3: "sp3"}
    for line in extract_lines:
        record, smiles, shifts = line.split("\t")
        if int(record) % 2:
            continue
        mol = Chem.MolFromSmiles(smiles)
        shifts = shifts.split(" ")
        if mol is None or mol.GetNumAtoms() != len(shifts):
            continue
        for atom in mol.GetAtoms():
            code = HYBRIDIZATION_CODES.get(atom.GetHybridization(), 0)
            shift = shifts[atom.GetIdx()]
            if (atom.GetSymbol() == "C" and atom.GetFormalCharge() == 0 and code
                    and shift not in (".", "-")):
                values[code].append(float(shift))
                molecules[code].add(int(record))
    lines = []
    for code in sorted(values):
        data = sorted(values[code])
        cuts = statistics.quantiles(data, n=1000, method="inclusive")
        numbers = [data[0]] + [cuts[c - 1] for c in HYBRIDIZATION_CUTS[:2]]
        numbers += [statistics.median(data)] + [cuts[c - 1] for c in HYBRIDIZATION_CUTS[2:]]
        numbers.append(data[-1])
        lines.append(f"hybridization-shift\t{names[code]}\t{len(data)}\t{len(molecules[code])}\t"
                     + "\t".join(f"{x:.1f}" for x in numbers))
    return lines


# ---- PubChem + RDKit oracle -------------------------------------------------------------


def pubchem_rows(path):
    rows = []
    for line in path.read_text().splitlines()[1:]:
        cid, smiles = line.split(",", 1)
        rows.append((int(cid), smiles.strip('"')))
    return sorted(rows)


def net_charge(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def cleaning_verdict(smiles):
    """RDKit's view of one PubChem record, for the cleaning check: 'mixture' (several
    components), 'ion' (charges that do not cancel), 'kept' (one neutral molecule, which
    may be charge-separated: nitro, N-oxide, zwitterion), or 'unreadable'."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unreadable"
    if len(Chem.GetMolFrags(mol)) > 1:
        return "mixture"
    return "ion" if net_charge(mol) != 0 else "kept"


def inchi_classes(mol):
    """Per atom, its class in standard InChI's equivalence layer (AuxInfo /E:), numbered
    by first appearance; atoms InChI lists in no class get one of their own. InChI applies
    charge and resonance normalisation (the two O of a nitro group are equivalent) and
    mobile-H (the two O of a carboxylic acid), so it may merge more than RDKit's ranks.
    None if InChI fails."""
    try:
        _, aux = Chem.MolToInchiAndAuxInfo(mol)
    except Exception:
        return None
    layers = dict(layer.split(":", 1) for layer in aux.split("/") if ":" in layer)
    if "N" not in layers:
        return None
    original = [int(x) - 1 for x in layers["N"].split(";")[0].split(",")]
    group = list(range(mol.GetNumAtoms()))
    for members in re.findall(r"\(([\d,]+)\)", layers.get("E", "").split(";")[0]):
        atoms = [original[int(x) - 1] for x in members.split(",")]
        for atom in atoms:
            group[atom] = atoms[0]
    first = {}
    return [first.setdefault(g, len(first)) for g in group]


def symmetry_classes(mol):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False))
    first = {}
    return [first.setdefault(rank, len(first)) for rank in ranks]


def arise_model_aromatic_atoms(mol):
    """Atoms of RDKit-aromatic rings that ARISE's model also calls aromatic: six-membered
    rings of C/N, and five-membered rings of four C plus one O, S or N."""
    atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() for a in ring_atoms):
            continue
        symbols = sorted(a.GetSymbol() for a in ring_atoms)
        six = len(ring) == 6 and set(symbols) <= {"C", "N"}
        five = len(ring) == 5 and symbols.count("C") == 4 and set(symbols) - {"C"} <= {"N", "O", "S"}
        if six or five:
            atoms.update(ring)
    return atoms


def aromaticity_differs(mol):
    aromatic = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    if aromatic - arise_model_aromatic_atoms(mol):
        return True
    for ring in mol.GetRingInfo().AtomRings():
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        for i in ring:
            for bond in mol.GetAtomWithIdx(i).GetBonds():
                other = bond.GetOtherAtomIdx(i)
                kind = bond.GetBondType()
                if kind == Chem.BondType.TRIPLE or (
                        kind == Chem.BondType.DOUBLE and other not in ring):
                    return True
    return False


def fragment_counts(mol):
    return [getattr(Fragments, name)(mol) for name in FRAGMENTS]


def fragment_patterns():
    """RDKit's published SMARTS of each fragment, read from its FragmentDescriptors.csv."""
    from rdkit import RDConfig
    patterns = {}
    with open(Path(RDConfig.RDDataDir) / "FragmentDescriptors.csv", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in FRAGMENTS:
                patterns[fields[0]] = Chem.MolFromSmarts(fields[2])
    if set(patterns) != set(FRAGMENTS):
        sys.exit(f"FragmentDescriptors.csv lacks {set(FRAGMENTS) - set(patterns)}")
    return patterns


FRAGMENT_PATTERNS = fragment_patterns()


# Which atoms of each RDKit fragment pattern carry the group's NMR hallmark (positions in
# the published SMARTS), for the atom-level check of chem::find_group's key atoms:
# fr_ketone [#6][CX3](=O)[#6] -> the C=O carbon is atom 1; fr_ester
# [#6][CX3](=O)[OX2H0][#6] -> the C=O carbon (1) and the O-C carbon (4); rings -> their
# carbons. A hallmark name per entry, as in chem::hallmarks.
KEY_ATOMS = {
    "fr_ketone": [("C=O carbon", [1])],
    "fr_aldehyde": [("C=O carbon", [0])],
    "fr_COO2": [("C=O carbon", [0])],
    "fr_ester": [("C=O carbon", [1]), ("O-C carbon", [4])],
    "fr_Al_OH": [("C-OH carbon", [0])],
    "fr_phenol": [("C-OH carbon", [1])],
    "fr_nitrile": [("C#N carbon", [1])],
    "fr_urea": [("C=O carbon", [0])],
    "fr_benzene": [("ring carbons", [0, 1, 2, 3, 4, 5])],
    "fr_pyridine": [("ring carbons", [1, 2, 3, 4, 5])],
    "fr_furan": [("ring carbons", [1, 2, 3, 4])],
    "fr_thiophene": [("ring carbons", [1, 2, 3, 4])],
    # fr_methoxy [OX2](-[#6])-[CH3]: the CH3 carbon is atom 2.
    "fr_methoxy": [("O-CH3 carbon", [2])],
}


def key_atoms(mol, patterns):
    """Per fragment (FRAGMENTS order) and hallmark: the sorted atoms RDKit's pattern
    matches at the KEY_ATOMS positions, as 'name=1,4' items joined by '/' ('-' if none)."""
    columns = []
    for name in FRAGMENTS:
        matches = mol.GetSubstructMatches(patterns[name], uniquify=True)
        items = []
        for hallmark, positions in KEY_ATOMS[name]:
            atoms = sorted({m[p] for m in matches for p in positions})
            items.append(hallmark.replace(" ", "_") + "=" + (",".join(map(str, atoms)) or "-"))
        columns.append("/".join(items))
    return ";".join(columns)


def read_smiles(smiles):
    """RDKit's reading of `smiles`, keeping explicit [H] atoms so the atom order matches
    ARISE's (ARISE keeps them as atoms and counts them toward their neighbour)."""
    params = Chem.SmilesParserParams()
    params.removeHs = False
    return Chem.MolFromSmiles(smiles, params)


# RDKit hybridization of a carbon as written in the oracle: 1 sp, 2 sp2, 3 sp3; 0 for
# any other atom or any other RDKit hybridization.
HYBRIDIZATION_CODES = {Chem.HybridizationType.SP: 1, Chem.HybridizationType.SP2: 2,
                       Chem.HybridizationType.SP3: 3}


def carbon_hybridizations(mol):
    return " ".join(str(HYBRIDIZATION_CODES.get(a.GetHybridization(), 0))
                    if a.GetSymbol() == "C" else "0" for a in mol.GetAtoms())


def oracle_line(smiles):
    mol = read_smiles(smiles)
    if mol is None or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    # Where RDKit's aromaticity model and ARISE's knowingly differ, the test skips the
    # comparisons that depend on it: rings ARISE does not model (azulene, azoles with two
    # heteroatoms, ...), and rings RDKit calls aromatic despite an exocyclic double bond
    # (quinoid tautomers) or a ring triple bond (arynes).
    differs = 1 if aromaticity_differs(mol) else 0
    hydrogens = " ".join(str(a.GetTotalNumHs(includeNeighbors=True)) for a in mol.GetAtoms())
    classes = " ".join(str(c) for c in symmetry_classes(mol))
    fragments = " ".join(str(c) for c in fragment_counts(mol))
    inchi = inchi_classes(mol)
    inchi_text = " ".join(str(c) for c in inchi) if inchi is not None else "-"
    keys = key_atoms(mol, FRAGMENT_PATTERNS)
    hybridizations = carbon_hybridizations(mol)
    return (f"{smiles}\t{hydrogens}\t{classes}\t{fragments}\t{differs}\t{inchi_text}\t{keys}"
            f"\t{hybridizations}")


def build_oracle(formulas):
    lines = []
    for formula in formulas:
        path = REPO / "tests" / "data" / f"pubchem_{formula}.csv"
        if not path.exists():
            download(PUBCHEM_URL.format(formula), path)
        seen = set()
        for _, smiles in pubchem_rows(path):
            if cleaning_verdict(smiles) != "kept" or smiles in seen:
                continue
            seen.add(smiles)
            line = oracle_line(smiles)
            if line is not None:
                lines.append(line)
        print(f"oracle {formula}: {len(seen)} distinct", file=sys.stderr)
    return lines


def build_pubchem_matches(sd_path, oracle):
    """Lines of nmrshiftdb2_pubchem_matches.tsv: every in-scope, InChI-checked nmrshiftdb2
    record that is one of the cleaned PubChem candidates of the oracle (same InChI skeleton),
    with that candidate's SMILES exactly as the oracle lists it. ARISE's tests rank all
    candidates of the record's formula against its assigned 13C spectrum (plan 3c)."""
    by_skeleton = {}
    for line in oracle:
        smiles = line.split("\t", 1)[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            by_skeleton.setdefault(inchi_skeleton(Chem.MolToInchi(mol)), smiles)
    rows = []
    for mol in Chem.ForwardSDMolSupplier(str(sd_path), removeHs=False, sanitize=True):
        if mol is None or not mol.HasProp("NMREDATA_ID") or not mol.HasProp("NMREDATA_INCHI"):
            continue
        match = re.search(r"DB_ID=(\d+)", mol.GetProp("NMREDATA_ID"))
        if not match or not mol.HasProp("NMREDATA_1D_13C"):
            continue
        listed = mol.GetProp("NMREDATA_INCHI")
        try:
            if not in_scope(mol) or not same_structure(mol, listed):
                continue
        except Exception:  # an RDKit failure on one record must not stop the build
            continue
        candidate = by_skeleton.get(inchi_skeleton(listed))
        if candidate is not None:
            rows.append((int(match.group(1)), candidate))
    rows = sorted(set(rows))  # a molecule with several spectra appears once
    print(f"nmrshiftdb2: {len(rows)} records are PubChem oracle candidates", file=sys.stderr)
    return [f"{db_id}\t{smiles}" for db_id, smiles in rows]


def build_cleaning_oracle(formulas):
    """RDKit's verdict on every record of the recorded PubChem formula searches."""
    lines = []
    for formula in formulas:
        for cid, smiles in pubchem_rows(REPO / "tests" / "data" / f"pubchem_{formula}.csv"):
            lines.append(f"{formula}\t{cid}\t{cleaning_verdict(smiles)}")
    return lines


# ---- benchmark case (plan 3d) -----------------------------------------------------------

# The benchmark cleans and prunes one real PubChem formula search of at least this many rows
# (ADR 0002's target; the chemist's decision, 2026-11-22: one formula, 10,000+ rows).
BENCHMARK_MIN_ROWS = 10000
BENCHMARK_MAX_FORMULAS = 25  # formulas tried before the build gives up
FORMULA_CACHE = PUBCHEM_CACHE / "formula"


def pubchem_formula_csv(formula):
    """PubChem's CID,ConnectivitySMILES CSV for `formula` (cached). PubChem sometimes answers
    'Status: 500'; the request is retried, and the build stops if it keeps failing, so the
    choice of formula never depends on a transient error."""
    FORMULA_CACHE.mkdir(parents=True, exist_ok=True)
    path = FORMULA_CACHE / f"{formula}.csv"
    for attempt in range(4):
        if path.exists() and path.read_text().startswith('"CID","ConnectivitySMILES"'):
            return path
        if attempt:
            time.sleep(5 * attempt)
        try:
            download(PUBCHEM_URL.format(formula), path)
        except OSError as error:
            print(f"PubChem {formula}: {error}", file=sys.stderr)
        time.sleep(0.25)  # PubChem allows at most 5 requests per second
    sys.exit(f"PubChem kept failing for {formula}; rerun later")


def benchmark_case(extract_lines):
    """The benchmark's formula, PubChem list and spectrum. Formulas of fully assigned
    nmrshiftdb2 extract records are tried from the most records to the fewest (ties by
    formula); the first whose PubChem search has at least BENCHMARK_MIN_ROWS rows is taken.
    Its record is the lowest-ID fully assigned one that is a cleaned PubChem candidate of
    that search (same InChI skeleton). Returns the header lines (each tried formula and its
    row count), the formula, and the data lines: record id, the candidate's PubChem SMILES,
    then one 'peak' line per distinct shift with the DEPT types of its carbons (RDKit H
    counts)."""
    by_formula = collections.defaultdict(list)
    for line in extract_lines:
        record, smiles, shifts = line.split("\t")
        values = shifts.split(" ")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or "-" in values:
            continue
        by_formula[rdMolDescriptors.CalcMolFormula(mol)].append((int(record), mol, values))
    order = sorted(by_formula, key=lambda f: (-len(by_formula[f]), f))
    tried = []
    for formula in order[:BENCHMARK_MAX_FORMULAS]:
        path = pubchem_formula_csv(formula)
        rows = pubchem_rows(path)
        tried.append(f"# tried {formula}: {len(by_formula[formula])} records, {len(rows)} PubChem rows")
        if len(rows) < BENCHMARK_MIN_ROWS:
            continue
        candidates = {}
        for _, smiles in rows:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None and cleaning_verdict(smiles) == "kept":
                candidates.setdefault(inchi_skeleton(Chem.MolToInchi(mol)), smiles)
        for record, mol, values in sorted(by_formula[formula], key=lambda r: r[0]):
            truth = candidates.get(inchi_skeleton(Chem.MolToInchi(mol)))
            if truth is None:
                continue
            peaks = collections.defaultdict(set)
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "C":
                    peaks[float(values[atom.GetIdx()])].add(
                        ["C", "CH", "CH2", "CH3", "CH4"][atom.GetTotalNumHs()])
            lines = [f"record\t{record}", f"structure\t{truth}"]
            lines += [f"peak\t{shift:g}\t{' '.join(sorted(types))}"
                      for shift, types in sorted(peaks.items())]
            return tried, formula, path, lines
    sys.exit("no nmrshiftdb2 formula has a PubChem search of "
             f"{BENCHMARK_MIN_ROWS}+ rows among the first {BENCHMARK_MAX_FORMULAS}")


# Compounds named by the chemist's question (charge-separated neutral molecules PubChem
# writes with formal charges); CID and formula come from PubChem, never typed.
NAMED_COMPOUNDS = ["nitrobenzene", "pyridine N-oxide", "trimethylamine N-oxide", "nitroethane"]
PUBCHEM_NAME_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/"
                    "MolecularFormula,ConnectivitySMILES/JSON")


def pubchem_by_name(name):
    """PubChem's one compound for `name`: (CID, formula, ConnectivitySMILES)."""
    url = PUBCHEM_NAME_URL.format(urllib.parse.quote(name))
    with urllib.request.urlopen(url, timeout=60) as response:
        props = json.load(response)["PropertyTable"]["Properties"]
    if len(props) != 1:
        sys.exit(f"PubChem gives {len(props)} compounds for {name!r}")
    p = props[0]
    return p["CID"], p["MolecularFormula"], p["ConnectivitySMILES"]


def named_compounds():
    lines = []
    for name in NAMED_COMPOUNDS:
        cid, formula, smiles = pubchem_by_name(name)
        lines.append(f"{name}\t{cid}\t{formula}\t{smiles}")
    return lines


# Molecules of the hand-written tests, which read them by name (the chemist's rule,
# 2026-11-21: no SMILES typed by hand in any test). Each is (name, PubChem CID, atom order):
# the structure is PubChem's (ConnectivitySMILES of the CID, or its isomeric SMILES for the
# CIDs in ISOMERIC_CIDS); the name is PubChem's title (or IUPAC name), with "(atom order n)"
# where one compound is used with several atom numberings. The atom order, when given,
# renumbers PubChem's atoms (atom k of the test molecule is PubChem's atom order[k]) so
# tests written against an earlier numbering keep it; the orders were computed once by
# RDKit substructure matching, and the builder checks that RDKit writes the renumbered
# molecule in exactly that order and that it is PubChem's structure (same InChIKey).
TEST_COMPOUNDS = [
    ('Borane', 6331, None),
    ('Bromofluoroiodomethane', 11481826, [2, 0, 1, 3]),
    ('Methane', 297, None),
    ('Formic Acid', 284, None),
    ('alpha-Glycine', 5257127, None),
    ('Butane-1,4-diolate', 22603739, None),
    ('3-Azaniumylpropan-1-olate', 58420841, None),
    ('Toluene', 1140, None),
    ('P-Xylene', 7809, None),
    ('Quinone', 4650, None),
    ('Fluoranthene', 9154, None),
    ('Toluene (atom order 2)', 1140, None),
    ('Biphenyl', 7095, None),
    ('Quinoline', 7047, None),
    ('Indole', 798, None),
    ('Benzofuran', 9223, None),
    ('Naphthalene', 931, None),
    ('Benzene', 241, None),
    ('1,3-Cyclohexadiene', 11605, [2, 3, 4, 5, 0, 1]),
    ('Pyridine', 1049, None),
    ('Pyridine N-Oxide', 12753, None),
    ('Cyclohexene', 8079, [3, 4, 5, 0, 1, 2]),
    ('Pyrrole', 8027, None),
    ('Furan', 8029, None),
    ('Thiophene', 8030, None),
    ('Cyclopentadiene', 7612, None),
    ('Gamma-Butyrolactone', 7302, None),
    ('Cyclohexane', 8078, None),
    ('Pyrrolidine', 31268, None),
    ('Tetrahydrofuran', 8028, None),
    ('Tetrahydrothiophene', 1127, None),
    ('2,3-Dihydrofuran', 70934, None),
    ('1,4-Dioxane', 31275, None),
    ('Isobutylene', 8255, [2, 1, 0, 3]),
    ('Linalool', 6549, [8, 7, 5, 6, 9, 4, 3, 2, 1, 0, 10]),
    ('Styrene', 7501, None),
    ('2-Ethenylfuran', 73881, None),
    ('1,3-Butadiene', 7845, None),
    ('Vinyl alcohol', 11199, None),
    ('Propane, 2-(ethenylsulfinyl)-', 5323637, [6, 5, 3, 4, 1, 0, 2]),
    ('Formaldehyde', 712, None),
    ('methylborane', 5326078, [1, 0]),
    ('Acetonitrile', 6342, None),
    ('Acetone', 180, None),
    ('3-[(Dimethylamino)methylidene]pentane-2,4-dione', 567003, [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
    ('Acetophenone', 7410, None),
    ('3-Buten-2-one, 4,4-bis(dimethylamino)-', 550970, None),
    ('Acetylacetone', 31261, [0, 1, 2, 3, 4, 6, 5]),
    ('Acetamide', 178, None),
    ('Dimethylacetamide', 31374, None),
    ('N-Acetylpyrrole', 521947, None),
    ('1-(2H-1,2,3-Triazol-2-yl)ethanone', 15939271, None),
    ('Acetohydrazide', 14039, None),
    ('Acetic Acid', 176, None),
    ('Acetic Anhydride', 7918, [0, 1, 2, 3, 4, 6, 5]),
    ('Iodobenzene diacetate', 76724, [0, 1, 2, 3, 4, 11, 12, 14, 13, 5, 6, 7, 8, 9, 10]),
    ('propan-2-one O-acetyl oxime', 2795071, [6, 4, 5, 3, 2, 1, 0, 7]),
    ('Neopentane', 10041, None),
    ('3-Diazo-2,2,4,4-tetramethylpentane', 638177, None),
    ('3,3-Dimethyl-2-butanol', 10045, [3, 2, 4, 5, 1, 0, 6]),
    ('2,3-Dimethyl-2-butene', 11250, [0, 1, 5, 2, 3, 4]),
    ('2,3-Pentadiene, 2-methyl-', 137807, [4, 3, 5, 2, 1, 0]),
    ('Schembl4340083', 12820292, [6, 5, 7, 4, 3, 1, 0, 2]),
    ('2-Propanamine, N-methyl-', 138743, [0, 1, 4, 2, 3]),
    ('2-Propanone, hydrazone', 78937, [0, 1, 4, 2, 3]),
    ('Acetone, oxime', 67180, [0, 1, 4, 2, 3]),
    ('2-Propanone, o-methyloxime', 137909, [0, 1, 5, 2, 3, 4]),
    ('DL-valine', 1182, [0, 1, 2, 3, 7, 4, 5, 6]),
    ('4-Methyl-2-Pentanol', 7910, None),
    ('Isopropanol', 3776, None),
    ('tert-butyl N-[methyl(propan-2-yl)-lambda4-sulfanylidene]carbamate', 21591014, [0, 1, 2, 3, 12, 4, 5, 6, 7, 8, 9, 10, 11]),
    ('Thioacetamide', 2723949, [0, 1, 3, 2]),
    ('1-Phenylethanol', 7409, [0, 1, 8, 2, 3, 4, 5, 6, 7]),
    ('1,1-Dimethylcyclohexane', 11549, [0, 1, 7, 2, 3, 4, 5, 6]),
    ('M-Xylene', 7929, [0, 1, 2, 3, 7, 4, 5, 6]),
    ('P-Xylene (atom order 2)', 7809, [0, 1, 2, 3, 4, 7, 5, 6]),
    ('2,5-Dimethylfuran', 12266, [0, 1, 2, 3, 4, 6, 5]),
    ('Toluene (atom order 3)', 1140, None),
    ('O-Xylene', 7237, None),
    ('Propylene Oxide', 6378, None),
    ('2-Methyl-2-butene', 10553, None),
    ('Methylketene', 138641, None),
    ('Acetaldehyde', 177, None),
    ('Propane', 6334, None),
    ('Methyl propionate', 11124, None),
    ('Schembl6298161', 132494933, None),
    ('Ethylbenzene', 7500, None),
    ('Triethylamine', 8471, None),
    ('Ethanol', 702, None),
    ('Ethyl Acetate', 8857, None),
    ('Ethyl Acetate (atom order 2)', 8857, [0, 1, 2, 3, 5, 4]),
    ('2-(Diethoxymethyl)-1,1-diethoxy-3-methylbutane', 638345, [0, 1, 2, 3, 15, 16, 17, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]),
    ('Diethyl Ether', 3283, None),
    ('Nitroethane', 6587, None),
    ('Trimethyl-lambda3-iodane', 20026379, None),
    ('Methylamine', 6329, None),
    ('Nitromethane', 6375, [0, 1, 3, 2]),
    ('Trimethylamine', 1146, None),
    ('1,1,3,3-Tetramethylguanidine', 66460, None),
    ('Tetramethylurea', 12437, None),
    ('4-(Dimethylamino)but-3-en-2-one', 566952, [6, 5, 7, 4, 3, 1, 0, 2]),
    ('Dimethylformamide', 6228, None),
    ('N,N-Dimethylthioformamide', 69794, None),
    ('Dimepranol', 37511, [4, 3, 5, 2, 1, 6, 0]),
    ('N-Nitrosodimethylamine', 6124, None),
    ('N,N-Dimethylmethanesulphonamide', 70191, None),
    ('1-Methylpyrrole', 7304, None),
    ('Dimethylamine', 674, None),
    ('Methomyl', 4109, [7, 6, 4, 5, 3, 2, 1, 0, 8, 9]),
    ('Aldicarb', 2086, [9, 8, 6, 7, 5, 4, 3, 1, 0, 2, 10, 11]),
    ('Methylamino methylcarbamate', 56997960, None),
    ('N-Methylacetamide', 6582, [4, 3, 1, 0, 2]),
    ('Methylhydrazine', 6061, None),
    ('Methylhydroxylamine', 11647, None),
    ('Dimethyl Ether', 8254, None),
    ('Methyl 5-methylhexa-2,3-dienoate', 14785091, [9, 8, 6, 7, 5, 4, 3, 1, 0, 2]),
    ('Carbamic acid, dimethyl-, methyl ester', 123254, [6, 5, 3, 4, 1, 0, 2]),
    ('Methyl hydrogen carbonate', 78579, None),
    ('Dimethyl Carbonate', 12021, None),
    ('Methyl Carbamate', 11722, [0, 1, 2, 4, 3]),
    ('Anisole', 7519, None),
    ('Methyl formate', 7865, None),
    ('Methoxyacetone', 22172, [5, 4, 3, 1, 2, 0]),
    ('Methoxyacetone (atom order 2)', 22172, [5, 4, 3, 1, 0, 2]),
    ('Methylphosphine', 68973, None),
    ('[Methyl(propan-2-yl)phosphoryl]benzene', 12847228, [5, 3, 4, 1, 0, 2, 6, 7, 8, 9, 10, 11]),
    ('Trimethylphosphine', 68983, None),
    ('Methanethiol', 878, None),
    ('Dimethyl Sulfone', 6213, None),
    ('Methanesulfonamide', 72879, [0, 1, 4, 2, 3]),
    ('Methanesulfinamide', 11115980, [0, 1, 3, 2]),
    ("a'-Dimethyl thiophene", 18717507, [0, 1, 6, 2, 3, 4, 5]),
    ('Schembl4657715', 57091962, None),
    ('Trimethylamine oxide', 1145, None),
    ('Methoxide', 3826051, None),
    ('Silaethane', 70434, None),
    ('Tetramethylsilane', 6396, None),
    ('Dichloromethane', 6344, [1, 0, 2]),
    ('Indole (atom order 2)', 798, [8, 7, 6, 4, 5, 0, 1, 2, 3]),
    ('Carbamic Acid', 277, [2, 0, 1, 3]),
    ('DL-Serine', 617, [5, 1, 0, 6, 2, 3, 4]),
    ('Urea', 1176, [2, 0, 3, 1]),
    ('Aniline', 6115, [6, 3, 2, 1, 0, 5, 4]),
    ('Formamide', 713, [2, 0, 1]),
    ('Glycine', 750, [4, 0, 1, 2, 3]),
    ('Water', 962, None),
    ('2-Pyridone', 8871, [3, 2, 1, 0, 6, 5, 4]),
    ('4-Hydroxypyridine', 12290, [6, 5, 4, 3, 2, 1, 0]),
    ('Succinic anhydride', 7922, [3, 2, 1, 0, 5, 6, 4]),
    ('2-Piperidone', 12665, [5, 4, 6, 0, 1, 2, 3]),
    ('2-Pyrrolidone', 12025, [3, 2, 1, 0, 5, 4]),
    ('2-Imidazolidinone', 8453, [4, 3, 2, 1, 0, 5]),
    ('2-Oxazolidinone', 73949, [4, 3, 5, 0, 1, 2]),
    ('Ethylene carbonate', 7303, [4, 3, 2, 1, 0, 5]),
    ('N-Nitrosopiperidine', 7526, [7, 6, 3, 2, 1, 0, 5, 4]),
    ('4-Hydroxybenzoic Acid', 135, [8, 6, 7, 5, 4, 3, 2, 9, 1, 0]),
    ('Benzoic Acid', 243, [8, 6, 7, 3, 2, 1, 0, 5, 4]),
    ('2-Naphthol', 8663, [10, 5, 4, 3, 2, 1, 0, 9, 8, 7, 6]),
    ('Phenol', 996, [6, 3, 2, 1, 0, 5, 4]),
    ('3-Hydroxypyridine', 7971, [6, 2, 3, 4, 5, 0, 1]),
    ('Cyclohexanol', 7966, [6, 3, 2, 1, 0, 5, 4]),
    ('Hexitol', 453, [6, 5, 4, 7, 3, 8, 2, 9, 1, 10, 0, 11]),
    ('Ethylene Glycol', 174, [2, 1, 0, 3]),
    ('Cyclohexanone Oxime', 7517, [5, 4, 3, 2, 1, 0, 7, 6]),
    ('Methyl isocyanide', 11646, [2, 1, 0]),
    ('Ammonium', 223, None),
    ('Bicyclo(2.2.1)heptane', 9233, None),
    ('(-)-2-Butanol', 84682, [3, 2, 4, 1, 0]),
    ('(+)-2-Butanol', 444683, [3, 2, 4, 1, 0]),
    ('trans-1,2-Difluoroethylene', 5365501, [2, 1, 0, 3]),
    ('cis-1,2-Difluoroethylene', 5462921, [2, 1, 0, 3]),
    ('Propane (1-13C)', 12314535, [2, 1, 0]),
    ('Carbon-13C', 105026, None),
    ('Methylselenol', 440764, None),
    ('Ethane', 6324, None),
    ('Methanol', 887, None),
    ('Thioanisole', 7520, None),
    ('Silane, methoxytrimethyl-', 15771, None),
    ('N,O-Dimethylhydroxylamine', 14232, None),
    ('Silane, methoxy-', 137466, None),
    ('Hydrogen Cyanide', 768, None),
    ('Methyl Ethyl Ketone', 6569, None),
    ('Butyric Acid', 264, None),
    ('Propylene', 8252, None),
    ('Dimethylethylamine', 11723, None),
    ('Trimethylphosphine oxide', 69609, None),
    ('Methylammonium', 644041, None),
    ('Ammonium Chloride', 25517, None),
    ('Sodium Acetate', 517045, None),
    # Substitute (2026-11-21, for the chemist's approval) for a hand-made zwitterion with
    # no PubChem entry: a quaternary N+ with two methyls and a stereocentre elsewhere.
    ('Lauryl Hydroxysultaine', 114474, None),
]

# Stereo and isotope labels matter to these tests: their structure is PubChem's isomeric SMILES.
ISOMERIC_CIDS = {84682, 444683, 5365501, 5462921, 12314535, 105026}
PUBCHEM_CID_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/"
                   "MolecularFormula,ConnectivitySMILES,SMILES/JSON")


def pubchem_by_cid(cid):
    """PubChem's (formula, ConnectivitySMILES, isomeric SMILES) for `cid`, cached."""
    cache = PUBCHEM_CACHE / f"cid_{cid}.json"
    if not cache.exists():
        cache.parent.mkdir(parents=True, exist_ok=True)
        download(PUBCHEM_CID_URL.format(cid), cache)
    p = json.loads(cache.read_text())["PropertyTable"]["Properties"][0]
    return p["MolecularFormula"], p["ConnectivitySMILES"], p["SMILES"]


def renumbered_smiles(name, smiles, order):
    """`smiles` itself when `order` is None; otherwise written by RDKit with its atoms in
    `order` (Kekule form); stops the build if RDKit would write them in another order or
    the structure changed."""
    if order is None:
        return smiles  # PubChem's own text, in PubChem's atom order
    mol = Chem.MolFromSmiles(smiles)
    if sorted(order) != list(range(mol.GetNumAtoms())):
        sys.exit(f"{name}: atom order is not a permutation of {mol.GetNumAtoms()} atoms")
    mol = Chem.RenumberAtoms(mol, order)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    text = Chem.MolToSmiles(mol, canonical=False, kekuleSmiles=True)
    written = [int(i) for i in mol.GetProp("_smilesAtomOutputOrder").strip("[],").split(",")]
    if written != list(range(mol.GetNumAtoms())):
        sys.exit(f"{name}: RDKit does not keep the atom order")
    if Chem.MolToInchiKey(Chem.MolFromSmiles(text)) != Chem.MolToInchiKey(
            Chem.MolFromSmiles(smiles)):
        sys.exit(f"{name}: the written SMILES is not PubChem's structure")
    return text


def charge_separated_nitrosamine(mol):
    """The chemist's rule (2026-11-21): writes the one N-N=O group of `mol` charge-separated.
    The amino N (single bonds only) gets a double bond to the nitroso N and a +1 charge; the
    N=O becomes a single bond and the O gets a -1 charge. Atoms, their order and their
    hydrogens are unchanged. Returns None unless there is exactly one such group."""
    pattern = Chem.MolFromSmarts("[NX3;!$(N=*);+0]-[NX2;+0]=[OX1;+0]")
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) != 1:
        return None
    amino, nitroso, oxygen = matches[0]
    edited = Chem.RWMol(mol)
    for atom in edited.GetAtoms():
        atom.SetNoImplicit(True)  # keep every atom's hydrogens as they are
        atom.SetNumExplicitHs(mol.GetAtomWithIdx(atom.GetIdx()).GetTotalNumHs())
    edited.GetBondBetweenAtoms(amino, nitroso).SetBondType(Chem.BondType.DOUBLE)
    edited.GetBondBetweenAtoms(nitroso, oxygen).SetBondType(Chem.BondType.SINGLE)
    edited.GetAtomWithIdx(amino).SetFormalCharge(1)
    edited.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
    result = edited.GetMol()
    Chem.SanitizeMol(result)
    return result


def uncharged_nitro(mol):
    """The chemist's rule (2026-11-21), the mirror of charge_separated_nitrosamine: writes
    the one charge-separated nitro group [N+](=O)[O-] of `mol` with two N=O double bonds
    and no charges. Atoms, their order and their hydrogens are unchanged. Returns None
    unless there is exactly one such group."""
    pattern = Chem.MolFromSmarts("[N+;X3](=[OX1;+0])-[OX1-]")
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) != 1:
        return None
    nitrogen, _, oxygen = matches[0]
    edited = Chem.RWMol(mol)
    for atom in edited.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(mol.GetAtomWithIdx(atom.GetIdx()).GetTotalNumHs())
    edited.GetBondBetweenAtoms(nitrogen, oxygen).SetBondType(Chem.BondType.DOUBLE)
    edited.GetAtomWithIdx(nitrogen).SetFormalCharge(0)
    edited.GetAtomWithIdx(oxygen).SetFormalCharge(0)
    result = edited.GetMol()
    result.UpdatePropertyCache(strict=False)  # RDKit's own sanitizing would recharge it
    return result


# Molecules made by a named rule from a PubChem compound (the chemist's decisions,
# 2026-11-21): (name, CID, rule).
RULE_COMPOUNDS = [
    ("N-Nitrosodimethylamine (charge-separated)", 6124, charge_separated_nitrosamine),
    ("Nitromethane (uncharged nitro)", 6375, uncharged_nitro),
]


def rule_compounds():
    lines = []
    for name, cid, rule in RULE_COMPOUNDS:
        formula, connectivity, _ = pubchem_by_cid(cid)
        source = Chem.MolFromSmiles(connectivity)
        made = rule(source)
        if made is None:
            sys.exit(f"{name}: rule {rule.__name__} does not apply exactly once")
        if (rdMolDescriptors.CalcMolFormula(made) != formula
                or Chem.GetFormalCharge(made) != 0
                or [a.GetTotalNumHs() for a in made.GetAtoms()]
                != [a.GetTotalNumHs() for a in source.GetAtoms()]
                or [a.GetSymbol() for a in made.GetAtoms()]
                != [a.GetSymbol() for a in source.GetAtoms()]):
            sys.exit(f"{name}: rule {rule.__name__} changed the formula, charge or atoms")
        if Chem.MolToInchiKey(made) != Chem.MolToInchiKey(source):
            sys.exit(f"{name}: InChI does not treat the rule's result as CID {cid}")
        Chem.Kekulize(made, clearAromaticFlags=True)
        text = Chem.MolToSmiles(made, canonical=False, kekuleSmiles=True)
        written = [int(i) for i in made.GetProp("_smilesAtomOutputOrder").strip("[],").split(",")]
        if written != list(range(made.GetNumAtoms())):
            sys.exit(f"{name}: RDKit does not keep PubChem's atom order")
        lines.append(f"{name}\t{cid}\t{rule.__name__}\t{connectivity}\t{text}\t"
                     f"{carbon_hybridizations(read_smiles(text))}")
    return lines


# Strings the SMILES-reader tests need that are not PubChem structures: notation variants
# and deliberately broken or unusual input (the chemist's decision, 2026-11-21: a named,
# RDKit-checked table instead of literals in the tests). Each is (name, string, compound):
# where `compound` names a TEST_COMPOUNDS entry, the string is another way of writing that
# compound, and the builder checks that RDKit reads it as the same structure (InChIKey).
SYNTAX_CASES = [
    ("empty input", "", None),
    ("two-digit ring closures", "C%10CCCCC%10", "Cyclohexane"),
    ("two-digit ring numbers 12 and 21", "C%12CC%21CC%12CC%21", "Bicyclo(2.2.1)heptane"),
    ("ring bond order on the opening digit", "C=1CCCCC1", "Cyclohexene"),
    ("ring bond order on the closing digit", "C1CCCCC=1", "Cyclohexene"),
    ("ring bond order on both digits", "C=1CCCCC=1", "Cyclohexene"),
    ("ring bond orders in conflict", "C=1CCCCC#1", None),
    ("explicit hydrogen atoms on methane", "[H]C([H])([H])[H]", "Methane"),
    ("explicit hydrogen atoms on one ethane methyl", "[H]C([H])([H])C", "Ethane"),
    ("explicit hydrogen atoms on propane C1", "[H]C([H])(C)C", "Propane"),
    ("explicit hydrogen atoms on methanol C", "[H]C([H])([H])O", "Methanol"),
    ("explicit hydrogen atom on methanol O", "[H]OC", "Methanol"),
    ("aromatic lowercase benzene", "c1ccccc1", "Benzene"),
    ("aromatic lowercase atom after a Kekule ring", "C1=CC=CC=C1[nH]", None),
    ("aromatic bond symbol", "CC:C", None),
    ("aromatic bracket atom", "[cH4]", None),
    ("atom class", "[CH4:1]", "Methane"),
    ("sodium", "C[Na]", None),
    ("silicon outside brackets", "CSi", None),
    ("quadruple bond symbol", "C$C", None),
    ("tetrahedral stereo class", "C[C@TH1H](O)C", None),
    ("two components", "CC.O", None),
    ("unclosed branch", "CC(C", None),
    ("branch closed without opening", "CC)C", None),
    ("empty branch", "C()C", None),
    ("branch before any atom", "(C)C", None),
    ("ring not closed", "C1CC", None),
    ("ring closed on its own atom", "C11", None),
    ("percent with one digit", "C%1C", None),
    ("two ring closures between the same atoms", "C12CC12", None),
    ("ring closure doubling a chain bond", "C1C1", None),
    ("ring closure onto the branch root", "CC(C1C)1", None),
    ("bond at the end", "CC=", None),
    ("bond at the start", "=CC", None),
    ("bond into an empty branch", "C(=)C", None),
    ("unclosed bracket", "C[CH4", None),
    ("empty bracket", "C[]C", None),
    ("doubly charged CH2 written ++", "[CH2++]", None),
    ("doubly charged CH2 written +2", "[CH2+2]", None),
    ("cationic CH2", "[CH2+]", None),
    ("CH in brackets between methyls", "C[CH]C", None),
    ("bare carbon in brackets", "[C]", None),
    ("CH2 radical", "[CH2]", None),
    ("CH2 radical end of propane", "[CH2]CC", None),
    ("O-CH fragment one hydrogen short", "CC(CO)O[CH]", None),
    ("methoxy radical", "C[O]", None),
    ("sulfur with three bonds", "C[S](C)C", None),
    ("sulfur with three bonds and one hydrogen", "C[SH](C)C", None),
    ("carbon with five bonds", "CC(C)(C)(C)C", None),
    ("CH5", "[CH5]", None),
    ("CF5", "FC(F)(F)(F)F", None),
    ("neutral nitrogen with four single bonds", "CN(C)(C)C", None),
    ("neutral nitrogen with five single bonds", "C[N](C)(C)(C)C", None),
    ("iodine with two bonds in brackets", "C[I]C", None),
    ("aromatic lowercase ethylbenzene", "CCc1ccccc1", "Ethylbenzene"),
    ("aromatic lowercase o-xylene", "Cc1ccccc1C", "O-Xylene"),
    ("aromatic lowercase m-xylene", "Cc1cccc(C)c1", "M-Xylene"),
    ("aromatic bond between aromatic atoms", "c:c", None),
    ("explicit single bond", "CC-O", "Ethanol"),
    ("mixture with a charged component", "CC(C)O.[C-]#[O+]", None),
    ("CH2 dication and two oxide anions, written ++", "[CH2++]C([O-])[O-]", None),
    ("CH2 dication and two oxide anions, written +2", "[CH2+2]C([O-])[O-]", None),
    ("CH2 dication and one oxide anion", "[CH2+2]CC[O-]", None),
    ("oxide dianion written --", "[O--]C[CH2+]", None),
    ("oxide dianion written -2, two cations", "[O-2]C[CH2+][CH2+]C", None),
    ("aromatic lowercase naphthalene", "c1ccc2ccccc2c1", "Naphthalene"),
    ("wildcard atom", "*C", None),
]


def syntax_cases(compound_smiles):
    """Rows of smiles_syntax_cases.tsv: name, string, the compound it writes ('-' if none),
    RDKit's verdict (accepted / rejected)."""
    names = [name for name, _, _ in SYNTAX_CASES]
    if len(set(names)) != len(names):
        sys.exit("SYNTAX_CASES lists a name twice")
    lines = []
    for name, text, compound in SYNTAX_CASES:
        mol = read_smiles(text) if text else None
        if compound is not None:
            if mol is None:
                sys.exit(f"syntax case {name!r}: RDKit cannot read it as {compound}")
            expected = Chem.MolToInchiKey(Chem.MolFromSmiles(compound_smiles[compound]))
            if Chem.MolToInchiKey(Chem.RemoveHs(mol)) != expected:
                sys.exit(f"syntax case {name!r} is not {compound}")
        verdict = "accepted" if mol is not None else "rejected"
        lines.append(f"{name}\t{text}\t{compound or '-'}\t{verdict}")
    return lines


def test_compounds():
    names = [name for name, _, _ in TEST_COMPOUNDS]
    if len(set(names)) != len(names):
        sys.exit("TEST_COMPOUNDS lists a name twice")
    lines = []
    for name, cid, order in TEST_COMPOUNDS:
        formula, connectivity, isomeric = pubchem_by_cid(cid)
        source = isomeric if cid in ISOMERIC_CIDS else connectivity
        smiles = renumbered_smiles(name, source, order)
        mol = read_smiles(smiles)
        lines.append(f"{name}\t{cid}\t{formula}\t{source}\t{smiles}\t"
                     f"{carbon_hybridizations(mol)}")
    return lines


# ---- hand-written spec cases and valences ---------------------------------------------

def build_spec_oracle(smiles_list):
    """RDKit's view of each SMILES the hand-written tests read (test_compounds.tsv,
    rule_compounds.tsv, smiles_syntax_cases.tsv): an oracle line (same columns as the
    PubChem oracle), or 'rejected' where RDKit refuses it (unparsable, valence error,
    radical)."""
    lines = []
    for smiles in sorted(set(smiles_list)):
        if "\t" in smiles or "\n" in smiles:
            continue
        line = oracle_line(smiles) if smiles else None
        lines.append(line if line is not None else f"{smiles}\trejected")
    return lines


# Organic-subset elements of OpenSMILES (implicit hydrogens from normal valences) and the
# elements ARISE accepts only in brackets.
ORGANIC_SUBSET = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
BRACKET_ONLY = ["Si", "Se"]


def valence_probes():
    """Probe molecules generated by rule, not chosen: each organic-subset element with
    0-6 methyl neighbours (implicit H), and each element in brackets with 0-4 explicit H and
    0-6 methyls. Columns: SMILES, RDKit's verdict ('ok' for a closed-shell molecule it
    accepts, else 'rejected'), RDKit's hydrogens on the probe atom (atom 0)."""
    smiles = []
    for element in ORGANIC_SUBSET:
        smiles += [element + "(C)" * k for k in range(7)]
    for element in ORGANIC_SUBSET + BRACKET_ONLY:
        for h in range(5):
            atom = f"[{element}H{h}]" if h else f"[{element}]"
            smiles += [atom + "(C)" * k for k in range(7) if h + k <= 6]
    # Charged probe atoms, each made net-neutral by a counter-charged group on its first
    # methyl (CC[O-] for a cation, CC[N+](C)(C)C for an anion) so they test valence, not
    # charge.
    for element in ORGANIC_SUBSET + BRACKET_ONLY:
        for sign, counter in (("+", "(CC[O-])"), ("-", "(CC[N+](C)(C)C)")):
            for h in range(5):
                atom = f"[{element}H{h}{sign}]" if h else f"[{element}{sign}]"
                smiles += [atom + counter + "(C)" * (k - 1) for k in range(1, 7) if h + k <= 6]
    lines = []
    for s in smiles:
        mol = read_smiles(s)
        if mol is None or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            lines.append(f"{s}\trejected\t-")
        else:
            lines.append(f"{s}\tok\t{mol.GetAtomWithIdx(0).GetTotalNumHs(includeNeighbors=True)}")
    return lines


def unsaturation(mol):
    """Rings (cyclomatic number: bonds - atoms + 1; RDKit's ring count is a symmetrized
    set that overcounts cages) plus pi bonds, on the Kekule graph."""
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    pi = sum(int(b.GetBondTypeAsDouble()) - 1 for b in kekule.GetBonds())
    return kekule.GetNumBonds() - kekule.GetNumAtoms() + 1 + pi


def lowest_valence_only(mol):
    """Whether every atom has the lowest valence RDKit's periodic table lists for it. A
    formula cannot tell S=O or P=O from S or P at their lowest valence, so the formula's
    degrees of unsaturation equal rings + pi bonds only for such molecules."""
    table = Chem.GetPeriodicTable()
    return all(a.GetTotalValence() == min(v for v in table.GetValenceList(a.GetAtomicNum())
                                          if v >= 0)
               for a in mol.GetAtoms())


def dou_line(smiles):
    """Molecular formula (RDKit) and rings + pi bonds of one sample molecule; None when
    the comparison does not apply (see lowest_valence_only)."""
    from rdkit.Chem import rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not lowest_valence_only(mol):
        return None
    return f"{rdMolDescriptors.CalcMolFormula(mol)}\t{unsaturation(mol)}"


# ---- ChEBI ---------------------------------------------------------------------------------


def chebi_terms(obo_path):
    """(is_a parents, name) of every ChEBI term."""
    parents = collections.defaultdict(set)
    names = {}
    term = None
    with gzip.open(obo_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("[Term]"):
                term = None
            elif line.startswith("id: CHEBI:"):
                term = int(line[len("id: CHEBI:"):].strip())
            elif line.startswith("name: ") and term is not None:
                names[term] = line[len("name: "):].strip()
            elif line.startswith("is_a: CHEBI:") and term is not None:
                parents[term].add(int(line[len("is_a: CHEBI:"):].split()[0]))
    return parents, names


def chebi_classes_of(parents):
    """ChEBI id -> the CHEBI_CLASSES among its is_a ancestors."""
    wanted = set(CHEBI_CLASSES)
    memo = {}

    def ancestors(term):
        if term not in memo:
            memo[term] = set()  # guards against cycles
            found = set()
            for parent in parents.get(term, ()):
                found.add(parent)
                found |= ancestors(parent)
            memo[term] = found
        return memo[term]

    return lambda term: ancestors(term) & wanted


def ring_counts(mol):
    return [getattr(Fragments, name)(mol) if pattern is None
            else len(mol.GetSubstructMatches(pattern, uniquify=True))
            for name, pattern in CHEBI_RINGS]


def verify_ring_patterns():
    """Checks the 'pubchem' rows of known_differences.tsv, which tie each ring pattern that
    is not an RDKit fragment to a PubChem compound: PubChem must still give that InChI and
    CID for the InChIKey, RDKit must derive the same InChIKey from the InChI, and the pattern
    must match the whole compound (every heavy atom), so it encodes exactly that ring."""
    rows = {}
    for line in (OUT / "known_differences.tsv").read_text().splitlines():
        fields = line.split("\t")
        if not line.startswith("#") and fields[0] == "pubchem":
            rows[fields[2]] = fields
    patterns = {name: pattern for name, pattern in CHEBI_RINGS if pattern is not None}
    if set(rows) != set(patterns):
        sys.exit(f"ring patterns {sorted(patterns)} but pubchem rows for {sorted(rows)}")
    for name, (_, inchi, _, reason) in rows.items():
        key = re.search(r"InChIKey ([A-Z]{14}-[A-Z]{10}-[A-Z])", reason).group(1)
        cid = int(re.search(r"PubChem CID (\d+)", reason).group(1))
        with urllib.request.urlopen(PUBCHEM_INCHIKEY_URL.format(key), timeout=60) as response:
            props = json.load(response)["PropertyTable"]["Properties"]
        if [(p["CID"], p["InChI"], p["InChIKey"]) for p in props] != [(cid, inchi, key)]:
            sys.exit(f"{name}: PubChem gives {props} for {key}, the table {cid} {inchi}")
        mol = Chem.MolFromInchi(inchi)
        if mol is None or Chem.MolToInchiKey(mol) != key:
            sys.exit(f"{name}: RDKit does not derive {key} from {inchi}")
        matches = mol.GetSubstructMatches(patterns[name])
        if len(matches) != 1 or len(matches[0]) != mol.GetNumHeavyAtoms():
            sys.exit(f"{name} does not match all of PubChem CID {cid}: {matches}")
        print(f"{name}: PubChem CID {cid}, {inchi}, {key} verified", file=sys.stderr)


def chebi_record(mol):
    """(Kekule SMILES, RDKit ring counts) of a ChEBI structure in ARISE's scope, or None."""
    if any(a.GetAtomicNum() == 0 or a.GetSymbol() not in ALLOWED or a.GetNumRadicalElectrons()
           for a in mol.GetAtoms()) or net_charge(mol):
        return None  # generic (R-group) structures, ions, radicals
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumHeavyAtoms() > CHEBI_MAX_HEAVY_ATOMS:
        return None
    rings = ring_counts(mol)
    heavy = Chem.RemoveHs(mol)
    Chem.Kekulize(heavy, clearAromaticFlags=True)
    smiles = Chem.MolToSmiles(heavy, kekuleSmiles=True, isomericSmiles=False, canonical=True)
    return smiles, rings


def chebi_inchis(structures_path):
    """ChEBI compound id -> the standard InChIs ChEBI lists for its structures."""
    import csv
    csv.field_size_limit(1 << 30)  # molfile columns are long
    inchis = collections.defaultdict(set)
    with gzip.open(structures_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["standard_inchi"]:
                inchis[int(row["compound_id"])].add(row["standard_inchi"])
    return inchis


def build_chebi(sdf_path, obo_path, structures_path):
    """(data lines, '# class <id>: <name>' header lines)."""
    parents, names = chebi_terms(obo_path)
    listed_inchis = chebi_inchis(structures_path)
    classes_of = chebi_classes_of(parents)
    members = collections.defaultdict(list)  # class -> [ChEBI id]
    records = {}
    charged_ids = set()
    with gzip.open(sdf_path, "rb") as handle:
        for mol in Chem.ForwardSDMolSupplier(handle):
            if mol is None:
                continue
            chebi_id = int(mol.GetProp("ChEBI ID").split(":")[1])
            classes = classes_of(chebi_id)
            if not classes:
                continue
            try:
                record = chebi_record(mol)
                if record is not None and not any(
                        same_structure(mol, inchi) for inchi in listed_inchis[chebi_id]):
                    DROPPED["chebi"] += 1
                    record = None
            except Exception:  # an RDKit failure on one record must not stop the build
                record = None
            if record is None:
                continue
            records[chebi_id] = (record[0], sorted(classes), record[1])
            if any(a.GetFormalCharge() for a in mol.GetAtoms()):
                charged_ids.add(chebi_id)
            for cls in classes:
                members[cls].append(chebi_id)
    # Neutral molecules and charge-separated ones (zwitterions, N-oxides, nitro) are
    # sampled separately, so widening the scope to the latter leaves the neutral sample
    # unchanged.
    sample = set()
    for cls in CHEBI_CLASSES:
        for charged, per_class in ((False, CHEBI_PER_CLASS), (True, CHEBI_CHARGED_PER_CLASS)):
            ids = sorted(i for i in members[cls] if (i in charged_ids) == charged)
            step = max(1, len(ids) // per_class)
            sample.update(ids[::step][:per_class])
        print(f"chebi CHEBI:{cls}: {len(members[cls])} molecules", file=sys.stderr)
    print(f"chebi: {len(records)} in scope, sampled {len(sample)}", file=sys.stderr)
    lines = [f"{i}\t{records[i][0]}\t{' '.join(str(c) for c in records[i][1])}"
             f"\t{' '.join(str(c) for c in records[i][2])}" for i in sorted(sample)]
    return lines, [f"# class {cls}: {names[cls]}" for cls in CHEBI_CLASSES]


def chebi_release(obo_path):
    with gzip.open(obo_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("data-version:"):
                return line.split(":", 1)[1].strip()
    return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--nmrshiftdb", type=Path,
                        default=Path.home() / ".cache" / "arise" / "nmrshiftdb2.nmredata.sd")
    parser.add_argument("--chebi", type=Path, default=Path.home() / ".cache" / "arise" / "chebi",
                        help="directory caching the ChEBI SDF and ontology downloads")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if not args.nmrshiftdb.exists():
        args.nmrshiftdb.parent.mkdir(parents=True, exist_ok=True)
        download(NMRSHIFTDB_URL, args.nmrshiftdb)
    refresh_license(OUT / "NMRSHIFTDB2_LICENSE.txt")
    nmr_lines, nmr_all = build_nmrshiftdb(args.nmrshiftdb)
    (OUT / "nmrshiftdb2_13c_all.tsv").write_text(header(
        ["# Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made",
         "# available here under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt,",
         f"# {LICENSE_URL}).",
         "# Every in-scope record with an assigned 13C spectrum that passes the InChI check,",
         "# for measuring each group's hallmark shifts (tools/measure_group_shifts).",
         "# Columns: nmrshiftdb2 id, Kekule SMILES, per-atom assigned 13C shift in ppm",
         "# ('.' = not carbon, '-' = carbon not listed in the assignment)."] + nmr_all,
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}"]))
    (OUT / "nmrshiftdb2_13c_sample.tsv").write_text(header(
        ["# Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made",
         "# available here under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt,",
         f"# {LICENSE_URL}).",
         "# Columns: nmrshiftdb2 id, Kekule SMILES, per-atom 13C peak index",
         "# ('.' = not carbon, '-' = carbon not listed in the assignment).",
         f"# Dropped {DROPPED['nmrshiftdb2']} records whose molfile does not give the InChI",
         "# nmrshiftdb2 lists for it (formula, connection and hydrogen layers)."] + nmr_lines,
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}"]))
    (OUT / "nmrshiftdb2_exchangeable.tsv").write_text(header(
        ["# Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made",
         "# available here under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt,",
         f"# {LICENSE_URL}).",
         "# SECONDARY evidence for hsqc::classify_protons: the HSQC below is derived from each",
         "# record's assignments (one cross peak per 1H peak on carbon, at that carbon's",
         "# assigned 13C shift), and the tests derive the D2O-shake list the same way. Only the",
         "# 1H shifts and which atom each H is on are measured data.",
         "# Columns: nmrshiftdb2 id, ARISE solvent, Kekule SMILES, 1H peaks as shift:kind",
         "# (kind C = H on carbon, X = H on O, N, S or Se), derived HSQC as h:c pairs."]
        + build_exchangeable_cases(args.nmrshiftdb),
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}"]))
    (OUT / "measurements.tsv").write_text(header(
        ["# Measured from nmrshiftdb2 (made available under the nmrshiftdb2 Database License,",
         "# NMRSHIFTDB2_LICENSE.txt) over in-scope records that pass the InChI check.",
         "# shift-range: nucleus, lowest and highest listed shift (ppm).",
         "# cn-shift: class (CN_CLASSES in the builder), n, median, 5th and 95th percentile",
         "# (statistics.quantiles, inclusive method) of the C=N carbon shift (ppm).",
         "# hybridization-shift: sp/sp2/sp3 (RDKit, uncharged carbons) over the 13C extract's",
         "# lines with even IDs: n carbons, molecules, min, 0.1st, 1st percentile, median,",
         "# 99th, 99.9th percentile, max (ppm; statistics.quantiles n=1000, inclusive)."]
        + measure_nmrshiftdb(args.nmrshiftdb) + measure_hybridization_shifts(nmr_all),
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}"]))
    # What the hand-written tests read by name: PubChem compounds, rule-made compounds and
    # the named syntax cases.
    compound_lines = test_compounds()
    rule_lines = rule_compounds()
    compound_smiles = {line.split("\t")[0]: line.split("\t")[4]
                       for line in compound_lines + rule_lines}
    syntax_lines = syntax_cases(compound_smiles)
    (OUT / "spec_rdkit_oracle.tsv").write_text(header(
        ["# RDKit's reading of every SMILES the hand-written tests read by name",
         "# (test_compounds.tsv, rule_compounds.tsv, smiles_syntax_cases.tsv).",
         "# Columns as pubchem_rdkit_oracle.tsv, or SMILES and 'rejected'.",
         "# fragments: " + " ".join(FRAGMENTS)]
        + build_spec_oracle(list(compound_smiles.values())
                            + [line.split("\t")[1] for line in syntax_lines]),
        [PUBCHEM_CID_URL, "SYNTAX_CASES in the builder"]))
    (OUT / "smiles_syntax_cases.tsv").write_text(header(
        ["# SMILES-reader test strings that are not PubChem structures (SYNTAX_CASES in the",
         "# builder; the chemist's decision, 2026-11-21), read by name in the tests.",
         "# Columns: name, string, the test compound it writes ('-' if none; checked by",
         "# InChIKey), RDKit's verdict (accepted / rejected)."] + syntax_lines,
        ["SYNTAX_CASES in the builder"]))
    (OUT / "valence_probes.tsv").write_text(header(
        ["# Rule-generated probe molecules (valence_probes in the builder). Columns: SMILES,",
         "# RDKit's verdict (ok = closed-shell molecule it accepts), RDKit's hydrogens on",
         "# atom 0."] + valence_probes(), ["RDKit's valence model"]))

    # The first three bring aromatic regioisomers; the last four bring the charge-separated
    # neutral molecules PubChem writes with formal charges (nitro, N-oxide, zwitterion,
    # isocyanide).
    tried, bench_formula, bench_csv, bench_lines = benchmark_case(nmr_all)
    shutil.copyfile(bench_csv, REPO / "tests" / "data" / f"pubchem_{bench_formula}.csv")
    (OUT / "benchmark_spectrum.tsv").write_text(header(
        ["# Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made",
         "# available here under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt,",
         f"# {LICENSE_URL}).",
         "# The benchmark case (bench/, plan 3d): the first nmrshiftdb2 formula, by number of",
         f"# fully assigned records, whose PubChem search has {BENCHMARK_MIN_ROWS}+ rows",
         "# (tests/data/pubchem_<formula>.csv), and its lowest-ID record that is a cleaned",
         "# candidate. Lines: formula; record (nmrshiftdb2 id); structure (the candidate's",
         "# PubChem SMILES); peak (shift in ppm, DEPT types of its carbons from RDKit H counts).",
         ] + tried + [f"formula\t{bench_formula}"] + bench_lines,
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}", PUBCHEM_URL.format(bench_formula)]))

    formulas = ["C4H8O2", "C8H8O2", "C10H12O2", "C6H5NO2", "C5H5NO", "C3H9NO", "C2H5NO2"]
    oracle = build_oracle(formulas)
    (OUT / "pubchem_cleaning_oracle.tsv").write_text(header(
        ["# RDKit's verdict on every record of the recorded PubChem formula searches",
         "# (tests/data/pubchem_<formula>.csv): mixture, ion (net charge not zero), kept",
         "# (one net-neutral molecule) or unreadable. Columns: formula, CID, verdict."]
        + build_cleaning_oracle(formulas + [bench_formula]),
        [PUBCHEM_URL.format(f) for f in formulas + [bench_formula]]))
    (OUT / "named_compounds.tsv").write_text(header(
        ["# Charge-separated neutral compounds looked up by name on PubChem (NAMED_COMPOUNDS",
         "# in the builder). Columns: name, CID, formula, ConnectivitySMILES."]
        + named_compounds(), [PUBCHEM_NAME_URL]))
    (OUT / "test_compounds.tsv").write_text(header(
        ["# Molecules the hand-written tests read by name (TEST_COMPOUNDS in the builder; the",
         "# chemist's rule, 2026-11-21: no hand-typed SMILES). Columns: name, CID, formula and",
         "# SMILES (from PubChem), the SMILES the tests read (RDKit's writing of PubChem's",
         "# structure in the test's atom order), RDKit's hybridization of each carbon of it",
         "# (1 sp, 2 sp2, 3 sp3; 0 for other atoms)."]
        + compound_lines, [PUBCHEM_CID_URL]))
    (OUT / "rule_compounds.tsv").write_text(header(
        ["# Molecules made by a named rule (RULE_COMPOUNDS in the builder) from a PubChem",
         "# compound, which the hand-written tests read by name. Columns: name, source CID,",
         "# rule, the source's ConnectivitySMILES (PubChem), the rule's result as the tests",
         "# read it, RDKit's hybridization of each carbon of it (1 sp, 2 sp2, 3 sp3; 0 other)."]
        + rule_lines, [PUBCHEM_CID_URL]))
    (OUT / "pubchem_rdkit_oracle.tsv").write_text(header(
        ["# Columns: SMILES, per-atom total H, per-atom symmetry class (first appearance),",
         "# RDKit fragment counts in the order of the 'fragments:' line, 1 if RDKit's",
         "# aromaticity knowingly differs from ARISE's (rings ARISE does not model, quinoid",
         "# tautomers, arynes), standard InChI's equivalence classes (AuxInfo /E:; '-' if",
         "# InChI fails), the key atoms RDKit's fragment patterns match (KEY_ATOMS; per",
         "# fragment 'hallmark=atoms' items, fragments separated by ';').",
         "# fragments: " + " ".join(FRAGMENTS),
         "# Last column: RDKit's hybridization of each carbon (1 sp, 2 sp2, 3 sp3; 0 for other",
         "# atoms or other hybridizations)."] + oracle,
        [PUBCHEM_URL.format(f) for f in formulas]))
    (OUT / "nmrshiftdb2_pubchem_matches.tsv").write_text(header(
        ["# Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made",
         "# available here under the nmrshiftdb2 Database License (NMRSHIFTDB2_LICENSE.txt,",
         f"# {LICENSE_URL}).",
         "# nmrshiftdb2 records that are cleaned PubChem candidates of pubchem_rdkit_oracle.tsv",
         "# (same InChI formula, connection and hydrogen layers). Columns: nmrshiftdb2 id, the",
         "# candidate's SMILES as the oracle lists it."]
        + build_pubchem_matches(args.nmrshiftdb, oracle),
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}"]
        + [PUBCHEM_URL.format(f) for f in formulas]))

    verify_ring_patterns()
    args.chebi.mkdir(parents=True, exist_ok=True)
    sdf = args.chebi / "chebi_lite_3_stars.sdf.gz"
    obo = args.chebi / "chebi_lite.obo.gz"
    structures = args.chebi / "structures.tsv.gz"
    for url, path in [(CHEBI_SDF_URL, sdf), (CHEBI_OBO_URL, obo),
                      (CHEBI_STRUCTURES_URL, structures)]:
        if not path.exists():
            download(url, path)
    chebi, class_names = build_chebi(sdf, obo, structures)
    samples = [line.split("\t")[1] for line in chebi + nmr_lines]
    dou = sorted({line for line in map(dou_line, samples) if line is not None})
    (OUT / "dou_checks.tsv").write_text(header(
        ["# Molecular formula and rings + pi bonds (RDKit, Kekule graph) of every molecule in",
         "# chebi_classes.tsv and nmrshiftdb2_13c_sample.tsv in which every atom has its",
         "# lowest RDKit valence (a formula cannot see S=O, P=O). Columns: formula, DoU.",
         "# Derived from ChEBI (CC BY 4.0) and nmrshiftdb2 (nmrshiftdb2 Database License,",
         "# NMRSHIFTDB2_LICENSE.txt)."] + dou,
        ["chebi_classes.tsv", "nmrshiftdb2_13c_sample.tsv"]))
    (OUT / "chebi_classes.tsv").write_text(header(
        [f"# Contains data from ChEBI ({chebi_release(obo)}), https://www.ebi.ac.uk/chebi/,",
         "# licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).",
         "# Kekule SMILES are generated by RDKit from ChEBI's molfiles; stereo is dropped.",
         "# Columns: ChEBI id, Kekule SMILES, the ChEBI classes of",
         "# docs/chem_group_definitions.md among the molecule's is_a ancestors, RDKit's",
         "# aromatic ring counts in the order of the 'rings:' line.",
         f"# Dropped {DROPPED['chebi']} in-scope records whose molfile does not give the",
         "# standard InChI ChEBI lists for the compound (structures.tsv).",
         "# ChEBI's own name of each class (chem_group_definitions.md must agree):"]
        + class_names + ["# rings: " + " ".join(name for name, _ in CHEBI_RINGS)] + chebi,
        [f"{CHEBI_SDF_URL} sha256={sha256(sdf)}", f"{CHEBI_OBO_URL} sha256={sha256(obo)}",
         f"{CHEBI_STRUCTURES_URL} sha256={sha256(structures)}"]))
    print(f"wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
