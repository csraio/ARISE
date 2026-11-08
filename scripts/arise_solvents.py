"""Matching a solvent name from a data source to one of ARISE's solvents (dev-time only).

Shared by scripts/build_bmrb_reference.py and scripts/build_reference_sets.py, which load it
by path (they run with python -I, which leaves the script directory off the import path).

A name matches an ARISE solvent (data/solvent_peaks.tsv) when PubChem gives exactly one
compound for it whose InChI skeleton (formula, connections, hydrogens) is that of the ARISE
solvent's PubChem CID, and it is deuterated: by the name itself (PubChem's compound carries
deuterium, as "D2O" or "CDCl3" do), by the source's own "deuterated" labelling, or, for a
name in UNLABELLED_DEUTERATED given with no labelling at all, by that rule.

ARISE never imports this module.
"""

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
HEADERS = {"Application": "ARISE reference builder"}

# Solvent names that count as deuterated when the source gives no isotopic labelling: NMR
# samples "in DMSO" are DMSO-d6 (the chemist's decision, 2026-11-04). Any other unlabelled,
# undeuterated name is still left out.
UNLABELLED_DEUTERATED = {"DMSO"}


def fetch(url, path=None, delay=0.25):
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


def pubchem_inchis(kind, key, cache):
    """The InChIs PubChem gives for a name or CID (kind "name" or "cid"), cached in
    `cache` (a directory)."""
    safe = re.sub(r"[^A-Za-z0-9]+", "_", str(key))
    path = cache / f"{kind}_{safe}.json"
    url = f"{PUBCHEM}/{kind}/{urllib.parse.quote(str(key), safe='')}/property/InChI/JSON"
    try:
        body = fetch(url, path)
    except urllib.error.HTTPError as error:
        if error.code != 404:
            raise
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'{"PropertyTable": {"Properties": []}}')
        body = path.read_bytes()
    return [p["InChI"] for p in json.loads(body)["PropertyTable"]["Properties"]]


def inchi_skeleton(inchi):
    """Formula, connections and hydrogens: the layers before stereo and isotopes."""
    layers = inchi.strip().split("/")
    kept = layers[:2]
    for layer in layers[2:]:
        if layer[:1] not in ("c", "h"):
            break
        kept.append(layer)
    return "/".join(kept)


def deuterated(inchi):
    """Whether the InChI's isotopic layer (from "/i" on, including its "/hD2"-style
    sublayer for exchangeable H) names deuterium."""
    at = inchi.find("/i")
    return at >= 0 and "D" in inchi[at:]


def arise_solvents(solvents_path, cache):
    """{InChI skeleton: ARISE solvent name} from data/solvent_peaks.tsv's PubChem CIDs."""
    found = {}
    for line in solvents_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        name, cid = line.split("\t")[:2]
        inchis = pubchem_inchis("cid", cid, cache)
        if len(inchis) != 1 or not deuterated(inchis[0]):
            sys.exit(f"{name}: PubChem CID {cid} is not one deuterated compound")
        found[inchi_skeleton(inchis[0])] = name
    return found


def identify(name, label, solvents, cache):
    """(ARISE solvent name, admitted by the UNLABELLED_DEUTERATED rule, None), or
    (None, False, the reason it is not one). `label` is the source's isotopic labelling
    text, or None when the source gives none."""
    inchis = pubchem_inchis("name", name, cache) if name else []
    if len(inchis) != 1:
        return None, False, "solvent name not one PubChem compound"
    labelled = label is not None and "deuterated" in label.lower()
    by_rule = label is None and name in UNLABELLED_DEUTERATED
    if not (deuterated(inchis[0]) or labelled or by_rule):
        return None, False, "solvent not deuterated"
    solvent = solvents.get(inchi_skeleton(inchis[0]))
    if solvent is None:
        return None, False, "solvent not one of ARISE's"
    return solvent, by_rule and not (deuterated(inchis[0]) or labelled), None
