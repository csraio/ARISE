# Proposed revisions to the ARISE proposal (for review)

These are the changes `docs/ARISE_Project_Proposal.html` needs to match ADR 0002. Once
approved, they will be applied to the HTML, and the PDF and DOCX regenerated. Nothing
in the proposal has been changed yet.

## 1. Summary and positioning
- **Now:** "a scoring and re-ranking engine … searches a library or re-ranks supplied candidates".
- **Proposed:** an engine that **follows an expert's workflow**:
  1. formula → every known structure from PubChem;
  2. cheap 1D NMR pruning (DEPT counts, symmetry, confirmed functional groups);
  3. shift prediction with the 3D anisotropy correction, and 2D scoring, only on the
     survivors.
  4. Re-ranking external lists (Sherlock, NMRfilter, .smi) stays as a secondary input.

## 2. Stage 1: Candidate sourcing
- **Remove:** library mode (a mass-indexed local candidate library built from a
  downloaded database).
- **Replace with:** a PubChem formula search (PUG REST `fastformula`), run by the web
  client. Then cleaning: deduplicate connectivities, and drop mixtures, charged
  species and radicals. Each removal is counted and explained.
- **Keep:** a high-resolution mass → candidate formulas as an optional front end.
- **Evidence:** live timings, 1.3–4.1 s for C4H8O2 to C15H24O (ADR 0002).

## 3. Stage 2: 1D pruning (replaces the "spin-system checklist")
- **Evidence used:**
  - C/CH/CH2/CH3 counts from 13C + optional DEPT-135/DEPT-90 (implemented: `nmr1d`);
  - symmetry-unique carbons versus 13C peaks;
  - functional groups suggested from shifts with their evidence, and confirmed by the
    chemist.
- **Only confirmed evidence prunes.** Every elimination is reported with its reason.

## 4. Stages 3–5
- **HOSE baseline:** unchanged in intent (an owned predictor, not nmrdb.org/SPINUS).
  Its detailed design is pending the author's input.
- **Anisotropy and 2D scoring:** unchanged in intent. Whether RDKit is used for
  conformers is decided at that stage.

## 5. Architecture
- **Remove:** "offline, air-gapped", and the local candidate library file.
- **Two deliverables:**
  - **Engine:** portable C++ (native CLI + WebAssembly), with no network or file access.
  - **Web client:** static files, with JavaScript for the UI and the PubChem request.
    No backend. It can be self-hosted or hosted publicly.
- **Dependencies:** none beyond GoogleTest for the pruning stages. No RDKit there.

## 6. Speed targets
- **Replace:** "Stages 1–2: 500 candidates < 50 ms".
- **With:** formula → pruned list in about 2–5 s end to end (dominated by PubChem);
  engine cleaning + pruning of 10,000 candidates < 100 ms.

## 7. Limitations (add)
- Only compounds present in PubChem can be found until built-in constrained isomer
  generation exists (future work).
- Results depend on PubChem's contents on the query date, so the retrieved list is
  saved with every result.

## 8. Roadmap
- **Reorder:**
  1. formula + PubChem handling;
  2. the in-house molecule reader;
  3. pruning filters;
  4. web client + CLI;
  5. HOSE baseline, anisotropy and 2D scoring.
- **Move:** "true de novo mode (surge)" stays as future work, now motivated by
  "compounds not in PubChem".
