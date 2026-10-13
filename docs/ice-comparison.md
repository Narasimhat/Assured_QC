# Comparison with ICE

Assured QC can run Synthego ICE on the same traces and put its scores beside its own. ICE is a comparator only: it changes no call, rank, score or figure. The comparison is made on the command line; ICE is Python software, so the browser version does not run it and nothing is uploaded.

## Running it

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

A case seen in practice is a compound heterozygote whose two deletions differ by one base. ICE placed only one of the two alleles and assigned the rest of the signal to small insertions, so its R squared was about 0.6 and its score fell well below Assured QC's. Fitting deletion sizes of 40 to 120 bases directly to the downstream peak heights of such a trace (with the control read as the reference) recovered both alleles and explained the trace far better than the alleles ICE placed; ICE's proposal set apparently does not contain the second allele, which is inferred from its fitted alleles and not from its code. The low ICE score is then a limit of ICE's model, not fewer knockout alleles. Such clones remain a reserve in the recommendation because the overlapping traces are hard to read; the comparison does not change that.

## Limits

SpCas9 only, one to three guides of 17 to 23 bases, a single read per sample, a donor of 20 to 300 bases; other samples are listed as not run with the reason. The ICE numbers depend on the ICE version, recorded in `run.json`.

ICE is used under its research and non-profit licence (copyright Synthego Corporation; the notice is retained in the ICE source you download) and is not redistributed with Assured QC.
