# Functional-group definitions (approved by the chemist, 2026-10-31)

ARISE's confirmable functional groups. A chemist confirms or rules out groups from the
spectra, and confirmed groups **eliminate** candidates. So the groups are fine-grained
motifs that NMR can tell apart. If they were coarse, a correct "no ester" would wrongly
eliminate a carbonate, whose C=O resonates elsewhere. A web page can still offer families
built from these exact groups (e.g. "any carbonyl derivative").

**Validation.** Each group maps to one ChEBI class, which is expert-curated (EBI; CC BY 4.0;
IDs verified via EBI OLS on 2026-10-31). Tests are automatic:
- whenever ChEBI asserts a class for a molecule, ARISE must find the group;
- ChEBI's annotation is incomplete, so a missing assertion is not checked;
- RDKit's published fragment library is a second, independent check, where it has an
  equivalent;
- every disagreement is an ARISE bug (fixed test-first) or a row in "Known differences"
  below, approved here.

**NMR hallmarks.** Each group's characteristic shift range will be **measured from
nmrshiftdb2's assigned spectra** (e.g. the 5th–95th percentile of the key carbon's shift),
not written from memory. 3b's suggestion rules use those measured ranges.

**Decided so far (user, 2026-10-31):**
- five-membered heteroaromatics are aromatic;
- formaldehyde counts as an aldehyde, formates as esters;
- carbonates are not esters, ureas and carbamates are not amides; each is its own group.

## Groups

| # | Group | ARISE counts | Includes | Excludes | ChEBI class | Key atom for the NMR hallmark |
|---|---|---|---|---|---|---|
| 1 | Carboxylic acid | C(=O)OH | formic acid | carboxylates (charged) | CHEBI:33575 carboxylic acid | C=O carbon; OH proton |
| 2 | Aldehyde | C(=O)H, carbonyl C with no single bond to O/N/S | formaldehyde | formic acid, formates, formamides | CHEBI:17478 aldehyde | C=O carbon; CHO proton |
| 3 | Ketone | C(=O) single-bonded to exactly two carbons | quinones, enones, aryl ketones | — | CHEBI:17087 ketone | C=O carbon |
| 4 | Carboxylic ester | C(=O)–O–C (that C not a carbonyl C), acyl side H or C | lactones, formates | carbonates, carbamates, anhydrides | CHEBI:33308 carboxylic ester | C=O carbon; O–C carbon |
| 5 | Carbonate ester | O–C(=O)–O with both O bonded to carbon | cyclic carbonates | — | CHEBI:46722 carbonate ester | C=O carbon |
| 6 | Carboxylic anhydride | C(=O)–O–C(=O) | cyclic anhydrides | — | CHEBI:35873 carboxylic anhydride | C=O carbons |
| 7 | Carboxamide | C(=O)–N, acyl side H or C | lactams, formamides | ureas, carbamates | CHEBI:37622 carboxamide | C=O carbon |
| 8 | Carbamate ester | O–C(=O)–N, the O also bonded to carbon or to an oxime N | cyclic carbamates, O-carbamoyl oximes (aldicarb, methomyl; user decision 2026-10-31) | carbamic acids | CHEBI:23003 carbamate ester | C=O carbon |
| 9 | Urea | N–C(=O)–N | cyclic ureas | — | CHEBI:47857 ureas | C=O carbon |
| 10 | Alcohol | OH on a saturated (sp3) carbon | hemiacetal OH | enols, phenols, acids | CHEBI:30879 alcohol | C–OH carbon |
| 11 | Enol | OH on a C=C carbon | — | phenols | CHEBI:33823 enol | C–OH carbon |
| 12 | Phenol | OH on an all-carbon aromatic ring | naphthols | hydroxypyridines (→ 13) | CHEBI:33853 phenols | C–OH carbon |
| 13 | Pyridone / hydroxypyridine | 2-/4-pyridone C=O in the ring, or OH on a pyridine carbon | both tautomers | — | CHEBI:38183 pyridone | C=O / C–OH carbon |
| 14 | Ether | O bonded to two carbons, neither a carbonyl C, O not in an aromatic ring | acetals, epoxides, aryl ethers | ester/anhydride O, furan O | CHEBI:25698 ether | O–C carbons |
| 15 | Amine | N bonded only to C/H by single bonds, not to C=O/C=S/S(=O), not in an aromatic ring | anilines | amides, sulfonamides, hydrazines, hydroxylamines, pyrrole N | CHEBI:32952 amine | N–C carbons |
| 16 | Sulfonamide | S(=O)2–N | — | — | CHEBI:35358 sulfonamide | — |
| 17 | Hydrazine | N–N single bond (not in an aromatic ring) | hydrazides' N–N | — | CHEBI:24631 hydrazines | — |
| 18 | Hydroxylamine | N–O single bond, N with only single bonds | — | oximes | CHEBI:24709 hydroxylamines | — |
| 19 | Imine | C=N, N bonded to no O or N | — | oximes, hydrazones | CHEBI:24783 imine | C=N carbon |
| 20 | Oxime | C=N–OH | — | — | CHEBI:25750 oxime | C=N carbon |
| 21 | Nitrile | C≡N | — | — | CHEBI:18379 nitrile | C≡N carbon |
| 22 | Alkene | C=C outside aromatic rings | enones, enols' C=C | aromatic ring bonds | CHEBI:78840 olefinic compound | C=C carbons |
| 23 | Benzene ring | aromatic all-carbon six-membered ring | each ring of naphthalene | pyridine | CHEBI:22712 benzenes | ring CH / C |
| 24 | Pyridine ring | aromatic six-membered ring with one N | — | — | CHEBI:26421 pyridines | ring carbons |
| 25 | Furan ring | aromatic five-membered ring with one O | — | — | CHEBI:24129 furans | ring carbons |
| 26 | Pyrrole ring | aromatic five-membered ring with one N | indoles | — | CHEBI:26455 pyrroles | ring carbons |
| 27 | Thiophene ring | aromatic five-membered ring with one S | — | — | CHEBI:26961 thiophenes | ring carbons |

## How the external checks map to the groups (user decisions, 2026-10-31)
- **ChEBI parent classes.** ChEBI files carbamate and carbonate esters under carboxylic
  ester, and pyridones under pyridines. When ChEBI asserts such a parent class, finding any
  of its ARISE groups passes.
- **ChEBI ring classes** (benzenes, pyridines, furans, pyrroles, thiophenes) are skeleton
  classes that include non-aromatic rings (butenolides, THF, dihydropyridines). They are
  checked only where RDKit independently sees the aromatic ring: its published fragments
  fr_benzene, fr_pyridine, fr_furan and fr_thiophene, and for pyrroles `n1cccc1`, written
  like RDKit's fr_furan `o1cccc1`. RDKit has no pyrrole fragment, so that pattern is tied in
  `known_differences.tsv` to the PubChem compound it encodes, 1H-pyrrole (PubChem CID 8027),
  by the InChI and InChIKey PubChem gives; the reference builder re-checks both against
  PubChem and RDKit and requires the pattern to match the whole compound.
- **RDKit fragments** are compared only where RDKit's published definition matches the
  group: fr_ketone, fr_aldehyde, fr_COO2 (acid), fr_ester, fr_Al_OH (alcohol + enol),
  fr_phenol, fr_nitrile, fr_urea, fr_benzene, fr_pyridine, fr_furan, fr_thiophene.
- **Oxime ethers** (C=N–O–C) get no group: their C=N carbon overlaps oximes' in
  nmrshiftdb2 (oximes: median 152.8 ppm, 5–95% 146.4–163.1, n = 98; oxime ethers: median
  156.4 ppm, 147.6–167.0, n = 160), so a separate group would not be told apart by shift.

## Known differences from the external checks
Every disagreement is listed by molecule in `tests/data/reference/known_differences.tsv`
(hand-written; each row names its category). A listed row that starts agreeing fails the
test, so the list cannot go stale. Categories, from the first full run (ChEBI release 256,
2,186 molecules; RDKit 2025.09.3, 4,655 PubChem candidates):

| Category | Rows | Source | Reason |
|---|---|---|---|
| formate | 51 | RDKit | fr_ester needs carbon on the acyl side; formates are esters (user decision) |
| chebi-error | 12 | ChEBI | ChEBI asserts a group whose atoms the structure lacks (e.g. "amine" with no N) |
| oxime-ether | 11 | ChEBI | ChEBI calls oxime ethers oximes and ethers; ARISE assigns no group (user decision) |
| tautomer | 10 | ChEBI | ChEBI classifies another tautomer (N-hydroxyguanidine, enol/keto forms, 4-iminopyridine) |
| urea-not-carboxamide | 5 | ChEBI | ureas are their own group (user decision) |
| pyridone-scope | 5 | ChEBI | ring with N=C, two ring C=O, or a dihydro ring: outside row 13 |
| carbonic-monoester | 4 | ChEBI | one O is OH; a carbonate ester needs carbon on both O |
| formamide | 3 | ChEBI | formamides are not aldehydes (user decision) |
| enol-as-alcohol | 2 | ChEBI | alcohol means OH on an sp3 carbon |
| enol-ester | 2 | ChEBI | enol means an OH on a C=C carbon |
| pyrrole-n, hydroxylamine-n, amide-n | 1 each | ChEBI | these N are not amines (row 15) |
| furan-o | 1 | ChEBI | an aromatic ring O is not an ether (user decision) |
| hydrazone | 1 | ChEBI | imines exclude C=N–N (row 19) |
| thiourea | 1 | ChEBI | urea needs C=O |
| ynol | 1 | RDKit | OH on an alkyne carbon: neither alcohol nor enol |
| bridged-benzene | 1 | RDKit | a bridged six-ring outside RDKit's smallest-ring set; ARISE finds the benzene |
| alpha-lactone | 1 | RDKit | fr_ester needs two distinct carbons around the ester O; ARISE counts it |

**Gaps these checks exposed** (no group today; add one if it matters): OH on a
heteroaromatic carbon other than pyridine (hydroxyfurans, pyrazolols), ynols, oxime ethers,
O-acyl oximes, hydrazones.

## Not covered yet
Sulfur groups beyond sulfonamide (thiol, sulfide, sulfoxide, sulfone), phosphorus groups,
halides (the formula already gives those), azoles with two heteroatoms (imidazole,
oxazole, thiazole, ...). Add when needed.
