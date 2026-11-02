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

**Charge-separated neutral molecules (user, 2026-11-01).** PubChem writes nitro groups,
N-oxides, zwitterions and isocyanides with formal charges that sum to zero; they are kept as
candidates. Groups count them by their neutral parent: an O⁻ on carbon counts as the OH it
came from (carboxylate as carboxylic acid, alkoxide as alcohol), an ammonium N as amine; an
N-oxide's N⁺–O⁻ is neither amine nor hydroxylamine. A carbonyl needs a neutral C=O (an
oxonium C=O⁺ or a carbonyl oxide C=O⁺–O⁻ is not one), a nitrile two neutral atoms, an imine
or oxime a neutral N.

## Groups

| # | Group | ARISE counts | Includes | Excludes | ChEBI class | Key atom for the NMR hallmark |
|---|---|---|---|---|---|---|
| 1 | Carboxylic acid | C(=O)OH, or C(=O)O⁻ in a net-neutral molecule (zwitterion) | formic acid, amino-acid zwitterions | carboxylate salts (ions are removed) | CHEBI:33575 carboxylic acid | C=O carbon; OH proton |
| 2 | Aldehyde | C(=O)H, carbonyl C with no single bond to O/N/S | formaldehyde | formic acid, formates, formamides | CHEBI:17478 aldehyde | C=O carbon; CHO proton |
| 3 | Ketone | C(=O) single-bonded to exactly two carbons | quinones, enones, aryl ketones | — | CHEBI:17087 ketone | C=O carbon |
| 4 | Carboxylic ester | C(=O)–O–C (that C not a carbonyl C), acyl side H or C | lactones, formates | carbonates, carbamates, anhydrides | CHEBI:33308 carboxylic ester | C=O carbon; O–C carbon |
| 5 | Carbonate ester | O–C(=O)–O with both O bonded to carbon | cyclic carbonates | — | CHEBI:46722 carbonate ester | C=O carbon |
| 6 | Carboxylic anhydride | C(=O)–O–C(=O) | cyclic anhydrides | — | CHEBI:35873 carboxylic anhydride | C=O carbons |
| 7 | Carboxamide | C(=O)–N, acyl side H or C | lactams, formamides | ureas, carbamates | CHEBI:37622 carboxamide | C=O carbon |
| 8 | Carbamate ester | O–C(=O)–N, the O also bonded to carbon or to an oxime N | cyclic carbamates, O-carbamoyl oximes (aldicarb, methomyl; user decision 2026-10-31) | carbamic acids | CHEBI:23003 carbamate ester | C=O carbon |
| 9 | Urea | N–C(=O)–N | cyclic ureas | — | CHEBI:47857 ureas | C=O carbon |
| 10 | Alcohol | OH (or O⁻ in a net-neutral molecule) on a saturated (sp3) carbon | hemiacetal OH | enols, phenols, acids | CHEBI:30879 alcohol | C–OH carbon |
| 11 | Enol | OH on a C=C carbon | — | phenols | CHEBI:33823 enol | C–OH carbon |
| 12 | Phenol | OH on an all-carbon aromatic ring | naphthols | hydroxypyridines (→ 13) | CHEBI:33853 phenols | C–OH carbon |
| 13 | Pyridone / hydroxypyridine | 2-/4-pyridone C=O in the ring, or OH on a pyridine carbon | both tautomers | — | CHEBI:38183 pyridone | C=O / C–OH carbon |
| 14 | Ether | O bonded to two carbons, neither a carbonyl C, O not in an aromatic ring | acetals, epoxides, aryl ethers | ester/anhydride O, furan O | CHEBI:25698 ether | O–C carbons |
| 15 | Amine | N bonded only to C/H by single bonds, not to C=O/C=S/S(=O), not in an aromatic ring | anilines | amides, sulfonamides, hydrazines, hydroxylamines, pyrrole N | CHEBI:32952 amine | N–C carbons |
| 16 | Sulfonamide | S(=O)2–N | — | — | CHEBI:35358 sulfonamide | — |
| 17 | Hydrazine | N–N single bond (not in an aromatic ring) | hydrazides' N–N | — | CHEBI:24631 hydrazines | — |
| 18 | Hydroxylamine | N–O single bond, N with only single bonds | — | oximes, N-oxides (N⁺–O⁻) | CHEBI:24709 hydroxylamines | — |
| 19 | Imine | C=N, N neutral and bonded to no O or N | — | oximes, hydrazones, iminium ions | CHEBI:24783 imine | C=N carbon |
| 20 | Oxime | C=N–OH | — | — | CHEBI:25750 oxime | C=N carbon |
| 21 | Nitrile | C≡N, both atoms neutral | — | isocyanides ([C⁻]≡[N⁺]) | CHEBI:18379 nitrile | C≡N carbon |
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
  nmrshiftdb2, so a separate group would not be told apart by shift. Measured by
  `scripts/build_reference_sets.py` (`measurements.tsv`; a test holds these rows to it):

  | Class | n | Median (ppm) | 5–95% (ppm) |
  |---|---|---|---|
  | oxime | 123 | 152.2 | 144.7–161.9 |
  | oxime-ether | 166 | 156.4 | 147.6–166.8 |

## Measured hallmark shifts (M3b; signed off by the chemist, 2026-11-01)
Measured over nmrshiftdb2's assigned 13C spectra with ARISE's own key atoms
(`chem::find_group`): records with even IDs; the odd half is kept back to check
suggestions. Each key atom of each record counts once. Generated by
`tools/measure_group_shifts.cpp` into `tests/data/reference/group_shifts.tsv`; a test
recomputes it from the committed extract and holds this table to it. Percentiles as
Python `statistics.quantiles(n=100)` (inclusive). Sulfonamide, hydrazine and
hydroxylamine have no carbon hallmark.

**Decisions (user, 2026-11-01):**
- A peak suggests a group when it falls in the group's **1st–99th percentile** window
  (about 98% of real occurrences). Suggestions never eliminate.
- The O–C (ester, ether) and N–C (amine) hallmarks are measured **separately for sp3 and
  unsaturated (aromatic, vinyl) carbons**; a peak in either window suggests the group.

**DEPT types (measured; awaiting the chemist's sign-off):** each row also lists the DEPT
types of the carbons it was measured on, read from their structures. A peak suggests the
group only if its possible DEPT types share one with the row's (a peak without DEPT
matches on shift alone), so DEPT separates, for example, an aldehyde CH from a ketone C
near 200 ppm. Because the types are measured, a group member missing from nmrshiftdb2 is
not covered: it has no formic acid (C=O is CH), methanol (C-OH is CH3), formaldehyde
(CH2) or HCN (CH).

**Hold-out check:** `tests/filter/holdout_test.cpp` runs the suggestions on the odd-ID
half, which the windows were not measured on. Per hallmark, the share of carbons whose own
peak falls in their window, and per group, the share of molecules with the group whose
spectrum suggests it, must reach 98% less three binomial standard errors
(sqrt(0.98 × 0.02 / n)); the test prints every rate.

**Using them (`filter`):** the engine holds no table of its own; the caller passes this
table's text to `filter::parse_hallmark_windows`, then `filter::suggest_groups` lists the
groups the 13C peaks suggest, each with its peaks and windows as evidence. Groups the
chemist confirms present or absent are applied by `filter::check_confirmed_groups`: a
candidate lacking a confirmed group, or having one ruled out, is Eliminated with one reason
per violated group. Nothing else eliminates.

| Group | Hallmark | DEPT types | n | Min | 1% | 5% | Median | 95% | 99% | Max |
|---|---|---|---|---|---|---|---|---|---|---|
| carboxylic acid | C=O carbon | {C} | 1334 | 154.2 | 160.4 | 163.1 | 172.7 | 182.7 | 185.5 | 189.2 |
| aldehyde | C=O carbon | {CH} | 566 | 92.4 | 175.9 | 178.3 | 191.0 | 202.7 | 205.1 | 207.0 |
| ketone | C=O carbon | {C} | 2981 | 115.8 | 169.0 | 174.9 | 195.8 | 213.4 | 218.8 | 225.8 |
| carboxylic ester | C=O carbon | {C, CH} | 3132 | 144.2 | 154.0 | 159.8 | 168.6 | 175.9 | 178.6 | 185.6 |
| carboxylic ester | O-C carbon (sp3) | {C, CH, CH2, CH3} | 2735 | 46.9 | 50.7 | 51.3 | 61.1 | 88.3 | 105.3 | 119.0 |
| carboxylic ester | O-C carbon (unsaturated) | {C, CH} | 390 | 77.7 | 128.3 | 137.5 | 150.8 | 161.0 | 174.9 | 193.0 |
| carbonate ester | C=O carbon | {C} | 34 | 146.2 | 146.6 | 148.2 | 153.6 | 155.8 | 156.6 | 156.9 |
| carboxylic anhydride | C=O carbon | {C} | 127 | 146.4 | 146.8 | 149.7 | 165.4 | 171.9 | 174.1 | 180.0 |
| carboxamide | C=O carbon | {C, CH} | 2067 | 140.8 | 152.2 | 155.8 | 167.6 | 177.5 | 181.6 | 192.0 |
| carbamate ester | C=O carbon | {C} | 273 | 146.4 | 147.7 | 149.1 | 155.0 | 160.2 | 162.4 | 165.0 |
| urea | C=O carbon | {C} | 298 | 144.3 | 146.9 | 149.1 | 152.8 | 164.3 | 168.5 | 175.9 |
| alcohol | C-OH carbon | {C, CH, CH2} | 3295 | 35.5 | 55.5 | 59.4 | 69.6 | 84.8 | 99.9 | 134.8 |
| enol | C-OH carbon | {C, CH} | 259 | 109.6 | 120.5 | 137.5 | 159.6 | 187.7 | 192.8 | 194.0 |
| phenol | C-OH carbon | {C} | 1884 | 123.4 | 133.0 | 140.8 | 154.1 | 164.7 | 166.7 | 169.3 |
| pyridone | C=O or C-OH carbon | {C} | 208 | 139.3 | 142.6 | 154.9 | 162.3 | 177.7 | 181.7 | 183.2 |
| ether | O-C carbon (sp3) | {C, CH, CH2, CH3} | 6977 | 43.6 | 48.7 | 53.8 | 66.4 | 103.3 | 115.7 | 135.8 |
| ether | O-C carbon (unsaturated) | {C, CH} | 4078 | 79.8 | 132.5 | 140.9 | 155.2 | 168.2 | 184.2 | 229.2 |
| amine | N-C carbon (sp3) | {C, CH, CH2, CH3} | 4142 | 23.3 | 30.0 | 36.2 | 51.1 | 69.9 | 90.2 | 127.0 |
| amine | N-C carbon (unsaturated) | {C, CH} | 2814 | 29.2 | 107.7 | 120.8 | 147.1 | 165.9 | 175.2 | 199.4 |
| imine | C=N carbon | {C, CH} | 1318 | 63.5 | 122.2 | 132.3 | 154.8 | 175.0 | 185.1 | 199.9 |
| oxime | C=N carbon | {C, CH} | 60 | 130.1 | 135.0 | 140.3 | 151.8 | 164.4 | 167.9 | 168.6 |
| nitrile | C#N carbon | {C} | 974 | 101.3 | 106.2 | 109.9 | 116.5 | 121.1 | 123.5 | 127.1 |
| alkene | C=C carbon | {C, CH, CH2} | 14223 | -3.0 | 76.6 | 97.9 | 131.1 | 163.2 | 178.9 | 218.7 |
| benzene ring | ring carbons | {C, CH} | 57787 | 65.1 | 99.8 | 110.3 | 128.7 | 154.3 | 162.3 | 176.5 |
| pyridine ring | ring carbons | {C, CH} | 5277 | 68.0 | 100.1 | 111.7 | 137.1 | 160.0 | 165.4 | 193.5 |
| furan ring | ring carbons | {C, CH} | 1755 | 61.2 | 96.1 | 104.1 | 134.2 | 158.9 | 163.7 | 167.0 |
| pyrrole ring | ring carbons | {C, CH} | 2569 | 57.8 | 93.3 | 102.5 | 125.4 | 141.5 | 150.0 | 158.0 |
| thiophene ring | ring carbons | {C, CH} | 1066 | 65.1 | 98.8 | 109.5 | 130.7 | 151.2 | 163.8 | 174.7 |

## Known differences from the external checks
Every disagreement is listed by molecule in `tests/data/reference/known_differences.tsv`
(hand-written; each row names its category). A listed row that starts agreeing fails the
test, so the list cannot go stale. Categories (ChEBI release 256, 2,321 molecules; RDKit
2025.09.3, 5,513 PubChem candidates; a test holds this table and these counts to the data
files):

| Category | Rows | Source | Reason |
|---|---|---|---|
| formate | 55 | RDKit | fr_ester needs carbon on the acyl side; formates are esters (user decision) |
| carbonyl-oxide | 6 | RDKit | C=[O+][O-] is a carbonyl oxide (Criegee intermediate), not a ketone or aldehyde; ARISE needs a neutral C=O |
| alkoxide-as-alcohol | 5 | RDKit | an O- on carbon counts as the OH of its neutral parent (user decision 2026-11-01); fr_Al_OH needs the H |
| carbamic-acid | 3 | RDKit | fr_COO2 counts the C(=O)OH or C(=O)O-; a carboxylic acid needs carbon or H on the acyl side, and carbamic acids are excluded from row 8 |
| oxonium-ylide | 2 | RDKit | an ylide/oxonium resonance artefact PubChem lists; ARISE needs a neutral C=O or OH where RDKit's pattern does not |
| alpha-lactone | 1 | RDKit | fr_ester needs two distinct carbons around the ester O; ARISE counts it |
| bridged-benzene | 1 | RDKit | a bridged six-ring outside RDKit's smallest-ring set; ARISE finds the benzene |
| ynol | 1 | RDKit | OH on an alkyne carbon: neither alcohol nor enol |
| chebi-error | 12 | ChEBI | ChEBI asserts a group whose atoms the structure lacks (e.g. "amine" with no N) |
| oxime-ether | 11 | ChEBI | ChEBI calls oxime ethers oximes and ethers; ARISE assigns no group (user decision) |
| tautomer | 10 | ChEBI | ChEBI classifies another tautomer (N-hydroxyguanidine, enol/keto forms, 4-iminopyridine) |
| pyridone-scope | 5 | ChEBI | ring with N=C, two ring C=O, or a dihydro ring: outside row 13 |
| urea-not-carboxamide | 5 | ChEBI | ureas are their own group (user decision) |
| carbonic-monoester | 4 | ChEBI | one O is OH; a carbonate ester needs carbon on both O |
| formamide | 3 | ChEBI | formamides are not aldehydes (user decision) |
| enol-ester | 2 | ChEBI | enol means an OH on a C=C carbon |
| amide-n | 1 | ChEBI | an amide N is not an amine (row 15) |
| enol-as-alcohol | 1 | ChEBI | alcohol means OH on an sp3 carbon |
| furan-o | 1 | ChEBI | an aromatic ring O is not an ether (user decision) |
| hydrazone | 1 | ChEBI | imines exclude C=N–N (row 19) |
| hydroxylamine-n | 1 | ChEBI | an N–O nitrogen is a hydroxylamine, not an amine (row 15) |
| pyrrole-n | 1 | ChEBI | the pyrrole N is aromatic, not an amine (row 15; user decision) |
| thiourea | 1 | ChEBI | urea needs C=O |
| ring-pattern | 1 | PubChem | RDKit has no pyrrole fragment; this pattern (written like fr_furan o1cccc1) is PubChem CID 8027 (1H-pyrrole), InChIKey KAESVJOAVNADME-UHFFFAOYSA-N |
| charged | 5 | chem tests | ARISE rejects charged species as candidates (molecule.hpp StructureIssue) |
| unsupported-notation | 5 | chem tests | ARISE reads Kekulé SMILES: no quadruple bonds, aromatic ':' bonds, lowercase aromatic atoms, @TH-style stereo classes or atom maps (molecule.hpp) |
| anhydride-not-ester | 2 | chem tests | fr_ester counts each C(=O)O of an anhydride; anhydrides are their own group (row 6) |
| alkoxide-as-alcohol | 2 | chem tests | an O- on carbon counts as the OH of its neutral parent (user decision 2026-11-01); fr_Al_OH needs the H |
| carbamic-acid | 1 | chem tests | fr_COO2 counts the C(=O)OH; a carboxylic acid needs carbon or H on the acyl side, and carbamic acids are excluded from row 8 |
| carbonic-monoester | 1 | chem tests | one O is OH; a carbonate ester needs carbon on both O |
| carbonyl-oxide | 1 | chem tests | C=[O+][O-] is a carbonyl oxide (Criegee intermediate), not a ketone or aldehyde; ARISE needs a neutral C=O |
| conflicting-ring-bond | 1 | chem tests | the two ends of a ring closure give different bond orders, an error in OpenSMILES; RDKit takes one of them |
| formaldehyde | 1 | chem tests | fr_aldehyde needs carbon on the carbonyl; formaldehyde is an aldehyde (user decision) |
| formate | 1 | chem tests | fr_ester needs carbon on the acyl side; formates are esters (user decision) |
| mixture | 1 | chem tests | ARISE reads one molecule; '.' is a ParseError (molecule.hpp) |
| unsupported-element | 1 | chem tests | outside ARISE's elements (molecule.hpp) |
| hypervalent-charged | 68 | valence probes | RDKit's valence lists for charged Cl, Br, I, P and Si (3-6 bonds on Cl+/Br+, ate complexes of I-, P-, Si-, P+ with 6 bonds) differ from ARISE's charge-shifted normal valences |
| small-sample | 1 | hold-out | measured on 127 carbons of 56 molecules, too few for its 1st-99th window to hold 98% of new data; cantharidin-type saturated cyclic anhydrides (175-178 ppm) fall above it (awaiting the chemist's sign-off) |

**Gaps these checks exposed** (no group today; add one if it matters): OH on a
heteroaromatic carbon other than pyridine (hydroxyfurans, pyrazolols), ynols, oxime ethers,
O-acyl oximes, hydrazones.

## Not covered yet
Sulfur groups beyond sulfonamide (thiol, sulfide, sulfoxide, sulfone), phosphorus groups,
halides (the formula already gives those), azoles with two heteroatoms (imidazole,
oxazole, thiazole, ...). Add when needed.
