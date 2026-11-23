# Functional-group definitions (approved by the chemist, 2026-10-31)

ARISE's confirmable functional groups. A chemist confirms or rules out groups from the
spectra, and confirmed groups **eliminate** candidates. So the groups are fine-grained
motifs. If they were coarse, a correct "no ester" would wrongly eliminate a carbonate,
which ARISE counts as its own group. A web page can still offer families
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
| 28 | Methoxy | CH3 single-bonded to an O that is bonded to exactly one other carbon; per OCH3 | aryl and aliphatic OMe, methyl esters, methyl carbonates and carbamates | methanol, O–CH3 on N, Si, S or P | CHEBI:25248 methyl ester; CHEBI:51683 methoxybenzenes | O–CH3 carbon |

**Methoxy (the chemist's decisions, 2026-11-22).** Added for 1H suggestions: a 3H singlet
of an OCH3 is one of the clearest 1H signs of a group. Methyl esters count as methoxy (their
OCH3 is the same signal); each OCH3 counts once, so dimethyl ether counts 2. ChEBI has no
methoxy class, so two narrower ChEBI classes check it one way: whenever ChEBI asserts
"methyl ester" (CHEBI:25248) or "methoxybenzenes" (CHEBI:51683), ARISE must find methoxy.
RDKit's fr_methoxy (`[OX2](-[#6])-[CH3]`, "Includes esters. Includes OMe groups on aromatic
rings.") checks it both ways; it counts dimethyl ether once (distinct atom sets), a
difference of counting only.

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
  fr_phenol, fr_nitrile, fr_urea, fr_benzene, fr_pyridine, fr_furan, fr_thiophene,
  fr_methoxy (methoxy; added 2026-11-22).
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
`tools/measure_group_shifts.cpp` into `data/group_shifts.tsv`; a test
recomputes it from the committed extract and holds this table to it. Percentiles as
Python `statistics.quantiles(n=100)` (inclusive). Sulfonamide, hydrazine and
hydroxylamine have no carbon hallmark.

**Decisions (user, 2026-11-01):**
- A peak suggests a group when it falls in the group's **1st–99th percentile** window
  (about 98% of real occurrences). Suggestions never eliminate.
- **Small groups (user, 2026-11-02):** a hallmark measured on fewer than 100 molecules
  (`filter::kMinMoleculesForPercentiles`) uses its **full measured range** (min–max)
  instead, because a percentile would rest on less than one molecule. Which rows this
  applies to follows from the Molecules column (today: carbonate ester, carboxylic
  anhydride and oxime).
- **ARISE ships this table (user, 2026-11-02)** and will be open source, as clause 4.5 of
  the nmrshiftdb2 licence requires of software relying on it. It lives in
  `data/group_shifts.tsv`, next to the licence, and is compiled into the engine
  (`filter::measured_hallmark_windows`).
- The O–C (ester, ether) and N–C (amine) hallmarks are measured **separately for sp3 and
  unsaturated (aromatic, vinyl) carbons**; a peak in either window suggests the group.

**DEPT types (measured; the chemist's choice, 2026-11-02):** each row also lists the DEPT
types of the carbons it was measured on, read from their structures. By default a peak
suggests the group only if its possible DEPT types share one with the row's (a peak
without DEPT matches on shift alone), so DEPT separates, for example, an aldehyde CH from a
ketone C near 200 ppm. Because the types are measured, a group member missing from
nmrshiftdb2 is not covered: it has no formic acid (C=O is CH), methanol (C-OH is CH3),
formaldehyde (CH2) or HCN (CH). The chemist can therefore switch DEPT off
(`filter::SuggestOptions::use_dept`, a tick-box in the web app), and every peak then
matches on shift alone.

**Hold-out check:** `tests/filter/holdout_test.cpp` runs the suggestions on the odd-ID
half, which the windows were not measured on. Per hallmark, the share of carbons whose own
peak falls in their window, and per group, the share of molecules with the group whose
spectrum suggests it, must reach 98% less three binomial standard errors
(sqrt(0.98 × 0.02 / n)); the test prints every rate.

**Using them (`filter`):** `filter::suggest_groups` lists the groups the 13C peaks
suggest, each with its peaks and windows as evidence.

**A lab's own windows:** `filter::lab_table_template` gives the built-in windows as a CSV
(columns `group, hallmark, dept, low_ppm, high_ppm`; only `group`, `low_ppm` and
`high_ppm` are required, names ignore case). The chemist edits it in a spreadsheet and
uploads it; `filter::parse_lab_windows` reads it and `filter::merge_windows` applies it. A
group the lab lists replaces all of ARISE's windows for that group; groups it leaves out
keep ARISE's. Groups the
chemist confirms present or absent are applied by `filter::check_confirmed_groups`: a
candidate lacking a confirmed group, or having one ruled out, is Eliminated with one reason
per violated group. Nothing else eliminates.

| Group | Hallmark | DEPT types | n | Molecules | Min | 1% | 5% | Median | 95% | 99% | Max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| carboxylic acid | C=O carbon | {C} | 1334 | 914 | 154.2 | 160.4 | 163.1 | 172.7 | 182.7 | 185.5 | 189.2 |
| aldehyde | C=O carbon | {CH} | 566 | 447 | 92.4 | 175.9 | 178.3 | 191.0 | 202.7 | 205.1 | 207.0 |
| ketone | C=O carbon | {C} | 2981 | 2169 | 115.8 | 169.0 | 174.9 | 195.8 | 213.4 | 218.8 | 225.8 |
| carboxylic ester | C=O carbon | {C, CH} | 3132 | 2416 | 144.2 | 154.0 | 159.8 | 168.6 | 175.9 | 178.6 | 185.6 |
| carboxylic ester | O-C carbon (sp3) | {C, CH, CH2, CH3} | 2735 | 2153 | 46.9 | 50.7 | 51.3 | 61.1 | 88.3 | 105.3 | 119.0 |
| carboxylic ester | O-C carbon (unsaturated) | {C, CH} | 390 | 314 | 77.7 | 128.3 | 137.5 | 150.8 | 161.0 | 174.9 | 193.0 |
| carbonate ester | C=O carbon | {C} | 34 | 32 | 146.2 | 146.6 | 148.2 | 153.6 | 155.8 | 156.6 | 156.9 |
| carboxylic anhydride | C=O carbon | {C} | 127 | 56 | 146.4 | 146.8 | 149.7 | 165.4 | 171.9 | 174.1 | 180.0 |
| carboxamide | C=O carbon | {C, CH} | 2067 | 1657 | 140.8 | 152.2 | 155.8 | 167.6 | 177.5 | 181.6 | 192.0 |
| carbamate ester | C=O carbon | {C} | 273 | 256 | 146.4 | 147.7 | 149.1 | 155.0 | 160.2 | 162.4 | 165.0 |
| urea | C=O carbon | {C} | 298 | 260 | 144.3 | 146.9 | 149.1 | 152.8 | 164.3 | 168.5 | 175.9 |
| alcohol | C-OH carbon | {C, CH, CH2} | 3295 | 2048 | 35.5 | 55.5 | 59.4 | 69.6 | 84.8 | 99.9 | 134.8 |
| enol | C-OH carbon | {C, CH} | 259 | 198 | 109.6 | 120.5 | 137.5 | 159.6 | 187.7 | 192.8 | 194.0 |
| phenol | C-OH carbon | {C} | 1884 | 971 | 123.4 | 133.0 | 140.8 | 154.1 | 164.7 | 166.7 | 169.3 |
| pyridone | C=O or C-OH carbon | {C} | 208 | 190 | 139.3 | 142.6 | 154.9 | 162.3 | 177.7 | 181.7 | 183.2 |
| ether | O-C carbon (sp3) | {C, CH, CH2, CH3} | 6977 | 3361 | 43.6 | 48.7 | 53.8 | 66.4 | 103.3 | 115.7 | 135.8 |
| ether | O-C carbon (unsaturated) | {C, CH} | 4078 | 2320 | 79.8 | 132.5 | 140.9 | 155.2 | 168.2 | 184.2 | 229.2 |
| amine | N-C carbon (sp3) | {C, CH, CH2, CH3} | 4142 | 1814 | 23.3 | 30.0 | 36.2 | 51.1 | 69.9 | 90.2 | 127.0 |
| amine | N-C carbon (unsaturated) | {C, CH} | 2814 | 1813 | 29.2 | 107.7 | 120.8 | 147.1 | 165.9 | 175.2 | 199.4 |
| imine | C=N carbon | {C, CH} | 1318 | 1137 | 63.5 | 122.2 | 132.3 | 154.8 | 175.0 | 185.1 | 199.9 |
| oxime | C=N carbon | {C, CH} | 60 | 57 | 130.1 | 135.0 | 140.3 | 151.8 | 164.4 | 167.9 | 168.6 |
| nitrile | C#N carbon | {C} | 974 | 811 | 101.3 | 106.2 | 109.9 | 116.5 | 121.1 | 123.5 | 127.1 |
| alkene | C=C carbon | {C, CH, CH2} | 14223 | 4706 | -3.0 | 76.6 | 97.9 | 131.1 | 163.2 | 178.9 | 218.7 |
| benzene ring | ring carbons | {C, CH} | 57787 | 7404 | 65.1 | 99.8 | 110.3 | 128.7 | 154.3 | 162.3 | 176.5 |
| pyridine ring | ring carbons | {C, CH} | 5277 | 912 | 68.0 | 100.1 | 111.7 | 137.1 | 160.0 | 165.4 | 193.5 |
| furan ring | ring carbons | {C, CH} | 1755 | 361 | 61.2 | 96.1 | 104.1 | 134.2 | 158.9 | 163.7 | 167.0 |
| pyrrole ring | ring carbons | {C, CH} | 2569 | 552 | 57.8 | 93.3 | 102.5 | 125.4 | 141.5 | 150.0 | 158.0 |
| thiophene ring | ring carbons | {C, CH} | 1066 | 235 | 65.1 | 98.8 | 109.5 | 130.7 | 151.2 | 163.8 | 174.7 |
| methoxy | O-CH3 carbon | {CH3} | 3423 | 2234 | 46.4 | 49.3 | 51.2 | 55.5 | 60.9 | 62.7 | 67.2 |

## Measured 1H windows (the 1H plan; the chemist's decisions, 2026-11-22)
Measured over nmrshiftdb2's assigned 1H spectra with ARISE's own key carbons
(`chem::find_group`), on records with even IDs: each key carbon's 1H shift is the mean of
the H on it (`tests/data/reference/nmrshiftdb2_1h_all.tsv`), split by the carbon's number of
H. Generated by `tools/measure_group_proton_shifts.cpp` into `data/group_proton_shifts.tsv`;
a test recomputes it from the committed extract and holds this table to it. The round-1
windows are also measured independently with RDKit's fragment patterns (`measurements.tsv`,
"proton-shift"), leaving out the records whose aromaticity the two models see differently.

**Decisions (the chemist, 2026-11-22):**
- Windows run from the 1st to the 99th percentile, or min–max from fewer than 100
  molecules, as for 13C. No windows for H on O, N or S (exchangeable protons).
- **Precision decides use.** One 1H peak in a crowded range says little, so each window also
  records its precision: of the even-ID molecules with an H of that count on any carbon
  inside the window, the share that have the group, against the base rate (the share of all
  molecules with the group).
- **Round 1 suggests only methoxy (O-CH3, CH3) and benzene ring (CH)**, the only windows
  right more often than not. The others are measured but not suggested; a later round enables
  any of them only if, on real PubChem searches, it keeps the true structure active and
  shrinks the active list, each signed off by the chemist.
- A group confirmed from 1H evidence will put contradicting candidates on the back burner,
  not eliminate them.

Columns: group, hallmark, H count, n carbons, molecules, min, 1st, 5th percentile, median,
95th, 99th percentile, max, window low, window high (ppm), hits, true hits, precision,
molecules with the group, molecules, base rate, round.

| Group | Hallmark | H | n | Molecules | Min | 1st | 5th | Median | 95th | 99th | Max | Low | High | Hits | True | Precision | With group | Molecules | Base | Round |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| carboxylic acid | C=O carbon | CH | 2 | 1 | 8.02 | 8.02 | 8.04 | 8.23 | 8.43 | 8.44 | 8.45 | 8.02 | 8.45 | 864 | 50 | 0.058 | 440 | 6408 | 0.069 | - |
| aldehyde | C=O carbon | CH | 180 | 168 | 5.05 | 6.30 | 8.23 | 9.86 | 10.53 | 11.05 | 11.47 | 6.30 | 11.05 | 3775 | 176 | 0.047 | 178 | 6408 | 0.028 | - |
| carboxylic ester | C=O carbon | CH | 7 | 7 | 8.03 | 8.03 | 8.03 | 8.05 | 8.13 | 8.13 | 8.13 | 8.03 | 8.13 | 346 | 59 | 0.171 | 1026 | 6408 | 0.160 | - |
| carboxylic ester | O-C carbon (sp3) | CH | 235 | 203 | 2.10 | 3.70 | 4.21 | 5.22 | 6.45 | 7.54 | 7.96 | 3.70 | 7.54 | 4579 | 740 | 0.162 | 1026 | 6408 | 0.160 | - |
| carboxylic ester | O-C carbon (sp3) | CH2 | 380 | 357 | 3.57 | 3.86 | 3.97 | 4.23 | 5.10 | 5.45 | 5.85 | 3.86 | 5.45 | 1645 | 499 | 0.303 | 1026 | 6408 | 0.160 | - |
| carboxylic ester | O-C carbon (sp3) | CH3 | 354 | 312 | 2.28 | 3.08 | 3.45 | 3.74 | 3.99 | 5.16 | 6.20 | 3.08 | 5.16 | 1121 | 361 | 0.322 | 1026 | 6408 | 0.160 | - |
| carboxylic ester | O-C carbon (unsaturated) | CH | 11 | 11 | 7.25 | 7.25 | 7.26 | 7.35 | 7.87 | 7.93 | 7.94 | 7.25 | 7.94 | 2333 | 301 | 0.129 | 1026 | 6408 | 0.160 | - |
| carboxamide | C=O carbon | CH | 20 | 20 | 6.59 | 6.83 | 7.78 | 8.26 | 8.83 | 9.21 | 9.30 | 6.59 | 9.30 | 3582 | 517 | 0.144 | 766 | 6408 | 0.120 | - |
| alcohol | C-OH carbon | CH | 578 | 395 | 2.33 | 3.04 | 3.35 | 4.08 | 5.54 | 6.25 | 6.58 | 3.04 | 6.25 | 2653 | 539 | 0.203 | 789 | 6408 | 0.123 | - |
| alcohol | C-OH carbon | CH2 | 367 | 316 | 3.06 | 3.27 | 3.41 | 3.75 | 4.75 | 5.15 | 5.35 | 3.27 | 5.15 | 2329 | 460 | 0.198 | 789 | 6408 | 0.123 | - |
| enol | C-OH carbon | CH | 3 | 3 | 6.50 | 6.53 | 6.67 | 8.18 | 8.83 | 8.89 | 8.90 | 6.50 | 8.90 | 3612 | 54 | 0.015 | 74 | 6408 | 0.012 | - |
| ether | O-C carbon (sp3) | CH | 924 | 529 | 2.33 | 2.82 | 3.34 | 4.59 | 6.23 | 6.68 | 7.81 | 2.82 | 6.68 | 3192 | 961 | 0.301 | 1418 | 6408 | 0.221 | - |
| ether | O-C carbon (sp3) | CH2 | 702 | 534 | 2.08 | 2.70 | 3.32 | 3.99 | 5.24 | 5.99 | 6.26 | 2.70 | 5.99 | 2904 | 837 | 0.288 | 1418 | 6408 | 0.221 | - |
| ether | O-C carbon (sp3) | CH3 | 718 | 513 | 2.56 | 3.02 | 3.24 | 3.79 | 4.01 | 4.19 | 5.47 | 3.02 | 4.19 | 1126 | 586 | 0.520 | 1418 | 6408 | 0.221 | - |
| ether | O-C carbon (unsaturated) | CH | 49 | 48 | 4.62 | 5.03 | 5.67 | 6.41 | 8.74 | 8.98 | 9.10 | 4.62 | 9.10 | 4705 | 1124 | 0.239 | 1418 | 6408 | 0.221 | - |
| amine | N-C carbon (sp3) | CH | 346 | 298 | 1.32 | 1.61 | 2.48 | 3.94 | 6.08 | 7.16 | 9.07 | 1.61 | 7.16 | 4260 | 841 | 0.197 | 1192 | 6408 | 0.186 | - |
| amine | N-C carbon (sp3) | CH2 | 628 | 426 | 1.35 | 2.05 | 2.36 | 3.10 | 4.61 | 5.28 | 6.26 | 2.05 | 5.28 | 3345 | 690 | 0.206 | 1192 | 6408 | 0.186 | - |
| amine | N-C carbon (sp3) | CH3 | 296 | 251 | 2.08 | 2.17 | 2.21 | 3.03 | 4.17 | 4.51 | 4.68 | 2.17 | 4.51 | 2036 | 476 | 0.234 | 1192 | 6408 | 0.186 | - |
| amine | N-C carbon (unsaturated) | CH | 222 | 197 | 4.73 | 5.31 | 6.30 | 8.09 | 9.54 | 9.78 | 10.05 | 5.31 | 9.78 | 4397 | 859 | 0.195 | 1192 | 6408 | 0.186 | - |
| imine | C=N carbon | CH | 157 | 147 | 6.26 | 6.88 | 7.53 | 8.36 | 9.53 | 10.11 | 11.22 | 6.88 | 10.11 | 3377 | 462 | 0.137 | 631 | 6408 | 0.098 | - |
| oxime | C=N carbon | CH | 4 | 4 | 7.29 | 7.29 | 7.30 | 7.66 | 8.41 | 8.47 | 8.49 | 7.29 | 8.49 | 2540 | 22 | 0.009 | 28 | 6408 | 0.004 | - |
| alkene | C=C carbon | CH | 3007 | 1622 | 0.88 | 4.30 | 5.15 | 6.66 | 8.63 | 9.54 | 11.22 | 4.30 | 9.54 | 4887 | 1841 | 0.377 | 1996 | 6408 | 0.311 | - |
| alkene | C=C carbon | CH2 | 335 | 311 | 3.85 | 4.05 | 4.65 | 5.20 | 6.24 | 6.55 | 7.70 | 4.05 | 6.55 | 1465 | 595 | 0.406 | 1996 | 6408 | 0.311 | - |
| benzene ring | ring carbons | CH | 9063 | 2461 | 4.35 | 6.30 | 6.63 | 7.38 | 8.32 | 8.90 | 10.52 | 6.30 | 8.90 | 3718 | 2756 | 0.741 | 3253 | 6408 | 0.508 | 1 |
| pyridine ring | ring carbons | CH | 679 | 241 | 5.58 | 6.33 | 6.94 | 8.12 | 9.23 | 9.77 | 10.22 | 6.33 | 9.77 | 3745 | 248 | 0.066 | 255 | 6408 | 0.040 | - |
| furan ring | ring carbons | CH | 201 | 85 | 4.90 | 5.22 | 5.87 | 7.04 | 8.04 | 8.22 | 8.70 | 4.90 | 8.70 | 4530 | 101 | 0.022 | 110 | 6408 | 0.017 | - |
| pyrrole ring | ring carbons | CH | 272 | 180 | 4.36 | 5.25 | 5.94 | 6.80 | 8.48 | 9.09 | 10.37 | 5.25 | 9.09 | 4395 | 248 | 0.056 | 257 | 6408 | 0.040 | - |
| thiophene ring | ring carbons | CH | 211 | 104 | 5.67 | 6.03 | 6.23 | 7.15 | 8.03 | 8.30 | 8.93 | 6.03 | 8.30 | 3790 | 114 | 0.030 | 133 | 6408 | 0.021 | - |
| methoxy | O-CH3 carbon | CH3 | 1082 | 807 | 2.28 | 3.04 | 3.28 | 3.77 | 4.00 | 4.31 | 6.20 | 3.04 | 4.31 | 1133 | 815 | 0.719 | 899 | 6408 | 0.140 | 1 |

## Known differences from the external checks
Every disagreement is listed by molecule in `tests/data/reference/known_differences.tsv`
(hand-written; each row names its category). A listed row that starts agreeing fails the
test, so the list cannot go stale. Categories (ChEBI release 256, 2,508 molecules; RDKit
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
| acetylated-phenol | 1 | ChEBI | ChEBI files a phenol's acetate ester under the parent phenol (syringaresinol diacetate); both ring OH are acetylated, so no phenol remains (row 12 needs an OH on the aromatic ring) |
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
| charged | 6 | chem tests | ARISE rejects charged species as candidates (molecule.hpp StructureIssue) |
| unsupported-notation | 10 | chem tests | ARISE reads Kekulé SMILES: no quadruple bonds, aromatic ':' bonds, lowercase aromatic atoms, @TH-style stereo classes, atom maps or '*' wildcard atoms (molecule.hpp) |
| anhydride-not-ester | 2 | chem tests | fr_ester counts each C(=O)O of an anhydride; anhydrides are their own group (row 6) |
| alkoxide-as-alcohol | 1 | chem tests | an O- on carbon counts as the OH of its neutral parent (user decision 2026-11-01); fr_Al_OH needs the H |
| carbamic-acid | 1 | chem tests | fr_COO2 counts the C(=O)OH; a carboxylic acid needs carbon or H on the acyl side, and carbamic acids are excluded from row 8 |
| carbonic-monoester | 1 | chem tests | one O is OH; a carbonate ester needs carbon on both O |
| carbonyl-oxide | 1 | chem tests | C=[O+][O-] is a carbonyl oxide (Criegee intermediate), not a ketone or aldehyde; ARISE needs a neutral C=O |
| methoxy-count | 1 | chem tests | fr_methoxy counts the O-CH3 atom sets of dimethyl ether once; ARISE counts each OCH3 (user decision 2026-11-22) |
| conflicting-ring-bond | 1 | chem tests | the two ends of a ring closure give different bond orders, an error in OpenSMILES; RDKit takes one of them |
| formaldehyde | 1 | chem tests | fr_aldehyde needs carbon on the carbonyl; formaldehyde is an aldehyde (user decision) |
| formate | 1 | chem tests | fr_ester needs carbon on the acyl side; formates are esters (user decision) |
| mixture | 4 | chem tests | ARISE reads one molecule; '.' is a ParseError (molecule.hpp) |
| unsupported-element | 1 | chem tests | outside ARISE's elements (molecule.hpp) |
| aromaticity | 1 | RDKit (1H) | the RDKit cross-check of the 1H windows leaves this nmrshiftdb2 record out: RDKit's model sees an aromatic six-ring of carbons in its fused ring system and ARISE's does not, and aromaticity_differs does not flag it |
| bmrb-assignment | 1 | BMRB | BMRB assigns 1H peaks that have HSQC cross peaks to H on O or N (guanosine, bmse000091): the entry's assignment contradicts its own HSQC list |
| hypervalent-charged | 68 | valence probes | RDKit's valence lists for charged Cl, Br, I, P and Si (3-6 bonds on Cl+/Br+, ate complexes of I-, P-, Si-, P+ with 6 bonds) differ from ARISE's charge-shifted normal valences |

**Gaps these checks exposed** (no group today; add one if it matters): OH on a
heteroaromatic carbon other than pyridine (hydroxyfurans, pyrazolols), ynols, oxime ethers,
O-acyl oximes, hydrazones.

## Not covered yet
Sulfur groups beyond sulfonamide (thiol, sulfide, sulfoxide, sulfone), phosphorus groups,
halides (the formula already gives those), azoles with two heteroatoms (imidazole,
oxazole, thiazole, ...). Add when needed.
