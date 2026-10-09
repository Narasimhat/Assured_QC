# Assured QC

Sanger-trace quality control for CRISPR edits in the browser: knockouts, deletions between cuts, SNP corrections, small tag knock-ins, and large tag or reporter knock-ins. The companion to the [Assured CRISPR Design](https://github.com/Narasimhat/Assured_CRISPR_Design) app: the design app writes what the correct allele should look like, this app checks whether your traces show it.

**Your files stay on your computer.** Traces are read and analysed inside the browser (in a web worker); nothing is uploaded and there is no server.

Status: 0.2.0, preview. See "Accuracy and limits" before using results to accept or reject clones.

## Workflows

| Workflow | Input | What you get |
| --- | --- | --- |
| Knockout (indels at one cut) | control + edited .ab1, guide(s) | editing %, knockout score, indel profile, clone genotype call |
| Deletion between two or three cuts | control + edited .ab1, guides | the above plus the deletion between the cuts |
| SNP correction (ssODN with silent blocking changes) | control + edited .ab1, design file or guides + donor (up to 300 nt) | intended edit separated from blocking-changes-only conversion, per-position readout of every designed change, clone genotype |
| Small tag knock-in (ssODN) | as for SNP | as for SNP |
| Large tag or reporter knock-in | junction or out-out reads + design file, optional band sizes | each read compared with the expected knock-in allele: junctions spanned, designed changes present, unexpected differences, out-out band classification |

To check the non-knock-in allele of a large knock-in clone, gel-purify the wild-type-sized band of the out-out PCR and run it as a Knockout sample.

Lab convention supported throughout: two guides per edit, with the donor carrying silent blocking changes for both guides.

## Use

**In the browser.** Drop a folder of .ab1 files (and, optionally, the design file exported by the design app, or a sample sheet) onto the start page. Files are paired automatically (control, edited, forward or reverse read) with a stated confidence for every pairing; the pairing table can be edited. The decision for each sample is accept, hold, re-sequence or reject, with the reasons, and the report, CSV, plate view, re-sequence list and a publication figure are one click each. The stepwise form is still available below the start page.

**From the command line** (same engine, parallel worker threads, no browser):

```
node cli/assured-qc.mjs run <folder> --design design.json --out results/      # or --sheet samples.csv, or --guides G1,G2 --donor <ssODN>
node cli/assured-qc.mjs figure --control control.ab1 --edited edited.ab1 --design design.json --out figure.svg
```

It writes `results.csv`, `results.json`, `report.html`, `plate.html`, `resequence.csv` and `run.json` (input hashes, parameters, engine version). 399 simulated pairs take 44 s on 11 threads (about 32,000 pairs per hour).

**Design file.** Export the design from the design app as a QC design file (JSON, schema `assured-qc-design/1`; see `docs/design-handoff.md`; needed for the reporter workflow), or type the guide sequence(s) and donor into the app. Other nucleases (SpCas9, SpCas9-NG, SaCas9, Cas12a, SpRY) are supported, and the cut can be inferred from the traces when no guide is given (simulation-validated only).

## What it does with a trace

The edited read is aligned to the design reference using the control read to remove the sequence-specific peak pattern. Within a window around the cut sites, the per-position base mixture of the edited trace is explained as a non-negative mixture of candidate alleles: wild type, indels of 1 to 30 bp at each cut, the deletion between cuts, the intended edit, and partial conversions (blocking changes without the full edit). The mixture weights are the allele fractions. Details in `docs/method.md`.

Unlike the headline editing number of other Sanger tools, the intended edit and the blocking-changes-only conversion are reported separately. A synthetic clone with half intended-edit alleles and half blocking-only alleles is reported here as about 50% intended edit and 47% blocking-only; ICE 1.2.0 reports 93% HDR for it (see `docs/validation.md`).

## Accuracy and limits

Held-out simulated test (traces built from real controls with known alleles; protocol and success criteria fixed before the run; `docs/scorecard.md`, `docs/results.md`): on 207 traces from loci not used in development, mean absolute error of the edited share is 1.9 points (95% CI 1.5-2.3) against 6.6 for ICE 1.2.0, 5.9 for peaksplit and 5.0 for Assured_TraceEdit; the intended-edit error is 4.2 against 7.0, 6.6 and 9.3, large deletions of 150 bp are recovered in 92% of traces (ICE 25%), and partial conversions are never reported as the intended edit (ICE and peaksplit always do). Criteria not met on this test: wild-type traces are called 2 or more points edited in 2 of 21 cases (the interval rule reports none of them as detected), the practical detection limit is about 5% (target 2%), HDR dilutions at 5-25% are 3 points off on average, and intended-edit intervals cover 81% of truths (target 90-98%; edited-share and knockout intervals cover 96% and 95%). The truth in this test is simulated, so these numbers are a bound on how well a tool can do on traces of this kind, not a measurement on your clones; the wet-lab mixture protocol (`docs/known-mixture-protocol.md`) has not been run. Only tools that can be scripted here were compared (ICE, peaksplit, Assured_TraceEdit); DECODR, TIDE, SeqScreener, EditR and Tracy were not.

On 172 laboratory pairs there is no sequencing-confirmed truth; the engine agrees with direct chromatogram readouts to 2.9 points for knockout clones (`docs/validation.md`). Earlier checks: 17 synthetic pairs used while the engine was built (1.5 points for editing %), and the ICE example traces (agreement within a few points, except a three-guide multiplex at 74% against ICE 63% and a 2-bp substitution at 34% HDR against ICE 44%).

- Detection limit is about 5% for an allele; Sanger traces cannot resolve smaller fractions. Every result carries an interval and a confidence tier (high, moderate, low) from the interval width; low-tier samples are not accepted automatically.
- A trace gives an unphased mixture per position. Alleles that differ only by combinations of changes, rather than by sequence, are reported as "also consistent with" where indistinguishable.
- The clone call assumes at most two alleles. A large deletion on one allele that removes a primer site shows up as apparent homozygosity; check with out-out PCR.
- Junction reads show whether a knock-in allele is correct, not how many copies or alleles are present.
- Trace quality limits everything. The app refuses swapped or mismatched control/edited pairs, reads that stop before the cut, and misaligned amplicons; it explains a control that is itself mixed (typically after a homopolymer run) rather than calling it the wrong design; it warns when the two traces already differ before the cut (discordance above 0.25, a provisional threshold from six real and several synthetic traces) and when the fit is poor.
- Real-trace accuracy beyond the ICE examples has to be confirmed with your own control/edited pairs of known outcome before this replaces your current check.

## Development

```
npm ci
npm test          # engine tests (node:test), no browser needed (107 tests)
npm run dev       # http://localhost:5174
npm run build
node scripts/validate.mjs out.csv   # own engine vs known mixtures vs ICE on the synthetic traces
python bench/simulate.py jobs_real.json out/ --seed N      # simulated truth from real controls
python bench/run_all.py jobs.json out/ --tools assured_qc,ice,peaksplit,traceedit
python bench/score.py jobs.json out/ prefix --split test    # scorer frozen in bench/FROZEN.json
```

The analysis modules in `src/` are plain ES modules with no browser or Node dependencies (`abif.js`, `align.js`, `decompose.js`, ...); `src/runJobs.js` is what the worker runs, `src/App.jsx` is the interface. `scripts/make_qc_fixtures.py` regenerates the synthetic traces (deterministic).

## Credits and licensing

The engine is written from scratch for this project. Synthego ICE (https://github.com/synthego-open/ice, Conant et al., CRISPR J 2022) was run only as a reference to compare results; none of its code is included, and its example traces are not redistributed. The repository licence has not been set yet.
