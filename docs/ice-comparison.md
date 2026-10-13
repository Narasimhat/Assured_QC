# Comparison with ICE

Assured QC puts the scores of ICE beside its own for the same samples. ICE is a comparator only: it changes no call, rank, score or figure. There are two ways to get ICE's numbers:

1. **Results from the ICE web tool** (browser and command line). You run the samples in the ICE web tool (ice.editco.bio), download the results and give the zip to Assured QC, which reads it and compares. Nothing is uploaded by Assured QC and ICE is not run again. This is the engine most people use today.
2. **ICE run by Assured QC** (command line only). The open-source ICE 1.2.0 from Synthego is run on the same traces. It is the older engine; see "Which ICE" below.

## Which ICE

The ICE version matters. The ICE web tool is now ICE 4.2, which EditCo describes as a new analysis engine; the open-source package is ICE 1.2.0. On a knockout plate of eight clones (two guides), Assured QC and ICE 4.2 agreed on all eight clones (KO score 100, R squared 0.94 to 1.00, the same deletion sizes in every clone), while ICE 1.2.0 scored two compound heterozygotes at 58 and 53 with R squared about 0.6. Compare with the web tool's results unless you need a local, scriptable run.

## Results from the ICE web tool

1. In the ICE web tool run the same control and edited traces with the same guides (for a point edit or tag, one run per donor, with the donor in the label, for example `A7_ssODN1`). Name each sample so that the label contains the sample name used in Assured QC (`cl4` is found in `cl4_MYD88_KO_1172`).
2. Download the results (a zip with `summary_<id>.csv` and one folder per sample).
3. Browser: under the results open "Compare with ICE" and add the zip (or the summary .csv); the zip can also be dropped with the other files. Command line: `node cli/assured-qc.mjs run traces --sheet sheet.csv --design design.html --out results --ice-results ICE_Results.zip` (the zip, its summary .csv, or the unzipped folder).

Only the summary is read: Label, ICE, KO-Score, KI-Score, R Squared, Guide Sequences, Notes and Indels (the allele sizes and shares). The layout was checked on an export of ICE 4.2; another layout is refused with a message, not guessed.

A label belongs to a sample when the sample's name occurs in it as whole words (`cl4` is not found in `cl40`; the longest name wins); a label that fits two samples equally well is not used, and labels that match no sample are listed. Samples without an ICE result are shown as "not run". The guides ICE was run with are compared with the design's, and a difference is written in the note, because the scores are then not like for like. ICE's own note, if any, is added.

The command line also writes `ice_results_used.json` (the rows read) and records the file name and SHA-256 of the export in `run.json`.

## Running ICE 1.2.0 on the command line

ICE is not bundled. Download the source of [synthego-open/ice](https://github.com/synthego-open/ice) (version 1.2.0 was used) and use a Python that has numpy, scipy, pandas and biopython.

    node cli/assured-qc.mjs run traces --sheet sheet.csv --design design_report.html --out results --ice /path/to/ice --ice-python /path/to/python

`--ice-python` defaults to `$AQC_ICE_PYTHON`, then `python` (`python3` outside Windows); `--ice-jobs N` sets the number of ICE processes (default: the number of analysis threads). ICE runs in separate processes because a process pool does not start on Windows.

Two steps, for a machine where Node cannot start Python, or to run ICE on a cluster node:

1. `run ... --ice-plan` writes `ice_jobs.json` (every ICE run: sample, donor, guides, trace paths).
2. `python scripts/run_ice_plan.py /path/to/ice results/ice_jobs.json results/ice_raw.json 6` runs ICE (anywhere Python and ICE are) and writes `ice_raw.json`.
3. `run ... --ice-raw results/ice_raw.json` reads it back. The results are accepted only for the same samples, guides, donors and traces (checked by SHA-256), so a report can also be remade later without running ICE again.

## What is written

- `ice_comparison.csv`: per sample the Assured QC score with its 95% interval, the ICE score, ICE's R squared, the difference, the status, the donor used, the alleles of both tools and the alleles ICE does not place.
- `ice_vs_assuredqc.svg` and `.pdf`: Assured QC against ICE, one point per clone; clones that agree are counted, the others are named.
- A section "Comparison with ICE" in `report.html`.
- `ice_raw.json`: what ICE returned, with the SHA-256 of every trace, and `ice` in `run.json` (ICE version, source, Python, thresholds, summary).

## How the two numbers are compared

- **Knockout and deletion:** the Assured QC KO score (frameshift or an indel of 21 bp or more, as in ICE) against ICE's KO score.
- **Point edit and small tag:** the Assured QC intended-edit share (its estimate, not the lower bound used for ranking) against ICE's HDR (knock-in) score. ICE takes one donor per run, so ICE is run once per donor and the best donor is shown, with the others listed in the note.
- **Agree:** within 15 points. **Differ:** 15 points or more with a sound ICE fit; look at the trace before relying on either.
- **ICE fit poor:** ICE's R squared below 0.8. ICE's own guidance is that such a fit should be read with caution, so its score is not comparable and is not counted as a disagreement. The note names the alleles Assured QC reads that ICE does not place.

A case seen in practice is a compound heterozygote whose two deletions differ by one base. ICE 1.2.0 placed only one of the two alleles and assigned the rest of the signal to small insertions, so its R squared was about 0.6 and its score fell well below Assured QC's. Fitting deletion sizes of 40 to 120 bases directly to the downstream peak heights of such a trace (with the control read as the reference) recovered both alleles and explained the trace far better than the alleles ICE 1.2.0 placed. The same traces in the ICE 4.2 web tool gave both alleles, a KO score of 100 and R squared 0.96 to 0.97, in agreement with Assured QC. The low ICE 1.2.0 score was therefore a limit of that engine, not fewer knockout alleles. Such clones remain a reserve in the recommendation because the overlapping traces are hard to read; the comparison does not change that.

## Limits

SpCas9 only, one to three guides of 17 to 23 bases, a single read per sample, a donor of 20 to 300 bases; other samples are listed as not run with the reason. The ICE numbers depend on the ICE version: it is recorded in `run.json` for ICE run by Assured QC, and for imported results the export itself is identified by its SHA-256. The import reads a zip with the standard library of the browser or of Node (no extra package), without zip64 or encryption.

The open-source ICE 1.2.0 is used under its research and non-profit licence (copyright Synthego Corporation; the notice is retained in the ICE source you download) and is not redistributed with Assured QC. The ICE web tool is used under its own terms by the person who runs it; Assured QC only reads the downloaded results.
