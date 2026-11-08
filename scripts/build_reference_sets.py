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
import sys
import urllib.parse
import urllib.request
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import Fragments

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
             "fr_nitrile", "fr_urea", "fr_benzene", "fr_pyridine", "fr_furan", "fr_thiophene"]

DEFINITIONS = REPO / "docs" / "chem_group_definitions.md"


def definition_classes():
    """The ChEBI class of each row of the signed-off groups table, in row order. Read from
    the doc so the IDs are typed in one place only (tests compare the C++ table with it)."""
    classes = []
    for line in DEFINITIONS.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\|\s*(\d+)\s*\|.*\|\s*CHEBI:(\d+) [^|]*\|", line)
        if match:
            if int(match.group(1)) != len(classes) + 1:
                sys.exit(f"{DEFINITIONS}: groups table row {match.group(1)} out of order")
            classes.append(int(match.group(2)))
    if not classes:
        sys.exit(f"{DEFINITIONS}: no groups table found")
    return classes


CHEBI_CLASSES = definition_classes()
# ChEBI's ring classes (furans, pyridines, ...) are skeleton classes that include
# non-aromatic rings (butenolides, dihydropyridines). ARISE's ring groups mean aromatic
# rings, so each ChEBI row also records RDKit's aromatic ring counts, and a ChEBI ring
# class is checked only where RDKit sees the ring (user decision 2026-10-31). Pyrroles have
# no RDKit fragment; "n1cccc1" is written like RDKit's fr_furan "o1cccc1", and
# known_differences.tsv ties it to the PubChem compound it encodes (verify_ring_patterns).
CHEBI_RINGS = [("fr_benzene", None), ("fr_pyridine", None), ("fr_furan", None),
               ("fr_thiophene", None), ("n1cccc1", Chem.MolFromSmarts("n1cccc1"))]
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
    return f"{smiles}\t{hydrogens}\t{classes}\t{fragments}\t{differs}\t{inchi_text}\t{keys}"


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


def build_cleaning_oracle(formulas):
    """RDKit's verdict on every record of the recorded PubChem formula searches."""
    lines = []
    for formula in formulas:
        for cid, smiles in pubchem_rows(REPO / "tests" / "data" / f"pubchem_{formula}.csv"):
            lines.append(f"{formula}\t{cid}\t{cleaning_verdict(smiles)}")
    return lines


# Compounds named by the chemist's question (charge-separated neutral molecules PubChem
# writes with formal charges); CID and formula come from PubChem, never typed.
NAMED_COMPOUNDS = ["nitrobenzene", "pyridine N-oxide", "trimethylamine N-oxide", "nitroethane"]
PUBCHEM_NAME_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/"
                    "MolecularFormula,ConnectivitySMILES/JSON")


def named_compounds():
    lines = []
    for name in NAMED_COMPOUNDS:
        url = PUBCHEM_NAME_URL.format(urllib.parse.quote(name))
        with urllib.request.urlopen(url, timeout=60) as response:
            props = json.load(response)["PropertyTable"]["Properties"]
        if len(props) != 1:
            sys.exit(f"PubChem gives {len(props)} compounds for {name!r}")
        p = props[0]
        lines.append(f"{name}\t{p['CID']}\t{p['MolecularFormula']}\t{p['ConnectivitySMILES']}")
    return lines


# ---- hand-written spec cases and valences ---------------------------------------------

# Test helpers in the hand-written tests (SPEC_DIRS) whose first argument is a SMILES literal.
SPEC_HELPERS = ("parse_smiles", "parse_ok", "parse", "expect_error", "groups", "environments",
                "carbon_class_sizes", "topological", "verdict", "exchangeable")
SPEC_DIRS = ("tests/chem", "tests/filter")


def spec_smiles():
    """Every SMILES literal the hand-written chem tests pass to a SPEC_HELPERS function,
    found by pattern so no list is typed by hand."""
    pattern = re.compile(r"\b(?:" + "|".join(SPEC_HELPERS) + r')\("((?:[^"\\]|\\.)*)"')
    found = set()
    for path in sorted(p for d in SPEC_DIRS for p in (REPO / d).glob("*_test.cpp")):
        for match in pattern.finditer(path.read_text(encoding="utf-8")):
            found.add(match.group(1).encode().decode("unicode_escape"))
    return sorted(found)


def build_spec_oracle():
    """RDKit's view of each spec SMILES: an oracle line (same columns as the PubChem
    oracle), or 'rejected' where RDKit refuses it (unparsable, valence error, radical)."""
    lines = []
    for smiles in spec_smiles():
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
    download(LICENSE_URL, OUT / "NMRSHIFTDB2_LICENSE.txt")
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
         "# (statistics.quantiles, inclusive method) of the C=N carbon shift (ppm)."]
        + measure_nmrshiftdb(args.nmrshiftdb),
        [f"{NMRSHIFTDB_URL} sha256={sha256(args.nmrshiftdb)}"]))
    (OUT / "spec_rdkit_oracle.tsv").write_text(header(
        ["# RDKit's reading of every SMILES literal in the hand-written tests",
         "# (tests/chem/*_test.cpp, tests/filter/*_test.cpp, found by pattern; SPEC_HELPERS",
         "# in the builder).",
         "# Columns as pubchem_rdkit_oracle.tsv, or SMILES and 'rejected'.",
         "# fragments: " + " ".join(FRAGMENTS)] + build_spec_oracle(),
        ["tests/chem/*_test.cpp", "tests/filter/*_test.cpp"]))
    (OUT / "valence_probes.tsv").write_text(header(
        ["# Rule-generated probe molecules (valence_probes in the builder). Columns: SMILES,",
         "# RDKit's verdict (ok = closed-shell molecule it accepts), RDKit's hydrogens on",
         "# atom 0."] + valence_probes(), ["RDKit's valence model"]))

    # The first three bring aromatic regioisomers; the last four bring the charge-separated
    # neutral molecules PubChem writes with formal charges (nitro, N-oxide, zwitterion,
    # isocyanide).
    formulas = ["C4H8O2", "C8H8O2", "C10H12O2", "C6H5NO2", "C5H5NO", "C3H9NO", "C2H5NO2"]
    oracle = build_oracle(formulas)
    (OUT / "pubchem_cleaning_oracle.tsv").write_text(header(
        ["# RDKit's verdict on every record of the recorded PubChem formula searches",
         "# (tests/data/pubchem_<formula>.csv): mixture, ion (net charge not zero), kept",
         "# (one net-neutral molecule) or unreadable. Columns: formula, CID, verdict."]
        + build_cleaning_oracle(formulas), [PUBCHEM_URL.format(f) for f in formulas]))
    (OUT / "named_compounds.tsv").write_text(header(
        ["# Charge-separated neutral compounds looked up by name on PubChem (NAMED_COMPOUNDS",
         "# in the builder). Columns: name, CID, formula, ConnectivitySMILES."]
        + named_compounds(), [PUBCHEM_NAME_URL]))
    (OUT / "pubchem_rdkit_oracle.tsv").write_text(header(
        ["# Columns: SMILES, per-atom total H, per-atom symmetry class (first appearance),",
         "# RDKit fragment counts in the order of the 'fragments:' line, 1 if RDKit's",
         "# aromaticity knowingly differs from ARISE's (rings ARISE does not model, quinoid",
         "# tautomers, arynes), standard InChI's equivalence classes (AuxInfo /E:; '-' if",
         "# InChI fails), the key atoms RDKit's fragment patterns match (KEY_ATOMS; per",
         "# fragment 'hallmark=atoms' items, fragments separated by ';').",
         "# fragments: " + " ".join(FRAGMENTS)] + oracle,
        [PUBCHEM_URL.format(f) for f in formulas]))

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
