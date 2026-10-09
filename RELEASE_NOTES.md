# Release notes

## 0.3.0 (preview)

- Engine: the read offset is re-anchored in the 20 bases before the cut when the two reads' base calls have drifted apart since the upstream alignment stretch (`src/offset.js`, option `refineOffset`). Fixes sizes that were one base out on a control read with extra base calls; never changes a pair whose offset already fits.
- Engine: unplanned base changes (a second base the design does not predict) are named with position and base (`src/unplanned.js`, option `unplanned`); results carry `unplanned`, a warning and a CSV column; an otherwise accepted clone is held (`holdOnUnplanned: false` to turn off).
- UI: the classic form's design box now takes an HTML design report (.html, .htm) as well as a .json design file (guides and recommended donor fill the form; before, only .json could be chosen). The Quick Start page already accepted reports.
- Bench: `bench/adapters_sangerdecon.py` (comparator adapter, separate from the frozen adapters) and `bench/PREREGISTRATION_0.3.0.md`.
- Evidence and limits: `docs/changes-0.3.0.md`. On 350 simulated traces (a fresh seed) the output is identical to 0.2.0.

## 0.2.0 (preview)

- Engine: unknown-base insertions to 10 bp; adaptive search for large deletions (31-250 bp, either cut or spanning), donor edit carrying an indel as its own allele; gated unpenalised refit and background baseline (only when the fit is good); joint fit of forward and reverse reads with a per-read agreement flag.
- Every result has profile-likelihood intervals (edited, wild type, intended edit, knockout score), a "could hide up to" bound, and a confidence tier derived from the interval width. The old any-warning flag did not predict errors (AUC 0.5); interval width does (AUC 0.8-0.9 on the development set).
- Decision layer: accept, hold, re-sequence or reject, with reasons; plate report (8x12) and re-sequence list. A held-out run found that an in-frame 18-base deletion between the cuts was accepted as a knockout; fixed and tested.
- Sample sheet and automatic pairing (167 of 172 real files correct, the 5 errors flagged low confidence), design-report import, zero-configuration start page with a worker pool, publication figure (SVG, PNG at 300 dpi).
- Command line (`cli/assured-qc.mjs`) with worker threads, run provenance (input SHA-256, parameters).
- Other nucleases (SpCas9, NG, SaCas9, Cas12a, SpRY), cut-site discovery from traces and a base-editing workflow (per-position conversion): checked on simulated traces only.
- Benchmark harness (`bench/`): simulator from real control peaks, adapters for ICE, peaksplit and Assured_TraceEdit, a scorer frozen before use, the known-mixture wet-lab protocol (not run). Held-out results in `docs/results.md`, including the criteria that were not met.
- 107 engine tests. The browser interface could not be built or viewed in the development environment (the bundler cannot run there); the production build and the tests pass in CI (Node 20), but nobody has looked at the interface in a browser yet.

## 0.1.0 (preview)

- Five workflows: knockout, deletion between cuts, SNP correction, small tag knock-in, large tag or reporter knock-in (junction reads and out-out band sizes).
- Design hand-off file (`assured-qc-design/1`) from the design app, or manual entry of guides and donor.
- Allele deconvolution with separate reporting of intended edit and blocking-changes-only conversion; clone genotype calls.
- HTML report and CSV table downloads; analysis runs in a web worker, files never leave the browser.
- Synthetic validation set and comparison against Synthego ICE 1.2.0 (`docs/validation.md`).
- Warning when the control and edited traces already differ before the cut (threshold 0.25, provisional); refusal of swapped, truncated and mismatched pairs. 68 engine tests.
- A control whose own signal is mixed (for example after a homopolymer run) is reported as such, with the position where the mixing starts and the run in the design, instead of as the wrong design (`src/diagnose.js`). Found on real reads from a project locus, where an (A)20 run upstream of the start codon turned every read mixed.
- Run on 172 laboratory control/edited pairs (12 sets, 160 clones, 12 pools) against ICE 1.2.0 and against direct chromatogram readouts (`docs/validation.md`, the clone-level table (not published)). Two failures found and fixed: a read with a poorly called start was refused although most of it aligned, and a complex pool was reported as 4% edited because the unexplained signal was renormalised away (new `unexplained` column in the table, report and CSV). 72 engine tests.
- A trace start that differs from the control by fewer than 10 bases is no longer listed as a warning (it was present in most pairs). `scripts/run_batch.mjs` runs a batch from the command line with the same engine and report as the web app.
- Known limits: three-guide multiplex headline differs from ICE on the one example tested; no sequencing-confirmed genotypes yet; poor traces remain the main source of error; at most two alleles per clone.
