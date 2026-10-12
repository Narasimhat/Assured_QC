# Scorecard (fixed before the comparison was run)

This page defines what "better" means for Assured QC against ICE, peaksplit and Assured_TraceEdit (Tracy has no Windows build; DECODR, TIDE/TIDER, SeqScreener and EditR are web-only and are not scripted). The scoring code is `bench/score.py`; its SHA-256 and the SHA-256 of the job file are in `bench/FROZEN.json`. Any later change to a definition is logged there as an amendment with its reason.

## Data
| Set | What | Truth | Use |
|---|---|---|---|
| D1 in-silico | 399 forward/reverse pairs from 14 real controls of 12 loci, 56 scenarios (`bench/simulate.py`, seed 20261009) | exact by construction | accuracy, detection limit, false positives, allele recovery |
| D2 real clones | 172 control/clone pairs, 12 loci | none (chromatogram readouts are checks, not truth) | regression guard: no engine change may worsen the readout agreement (`docs/validation.md`) |
| D3 known mixtures | wet-lab series in `docs/known-mixture-protocol.md` | bench molar fractions | confirmatory; not yet run |
| D4 L12 clones | 8 clones | junctions read base by base | allele identity |

Development / test split by locus, fixed now: **development** loci are L01, L02, L04 (knockout sets) and L07, L09 exon 20, L11 (SNP sets); **test** loci are L03, L05 exon 4, L06, L08, L10, L05 exon 25. Engine changes are tuned on development loci only. Headline numbers are reported on the test loci, on a fresh generator seed (20261010) that is not used during development, and again on all loci.

## Metrics (per tool; failures to analyse are reported, never dropped silently)
1. **Edited-fraction error.** Mean absolute error (MAE, percentage points) of the non-wild-type share, with a bootstrap 95 % interval over traces, overall and in truth bands 0, (0,5], (5,25], (25,75], (75,100]. Two versions: over analysed traces, and over all traces with a failure counted as "wild type called" (edited = 0).
2. **False positives.** Share of wild-type replicates with edited >= 2.0 points.
3. **Detection limit.** A trace is detected if edited >= 2.0 points. Sensitivity is computed per level in the graded series (two-cut deletion 2-100 %, HDR 2-100 %, +1 insertion 2/5/10 %). Limit = lowest level with sensitivity >= 90 % while the false-positive share is <= 5 %.
4. **Intended-edit error (SNP loci).** MAE of intended-edit share against truth, and the **false-intended rate**: traces with truth < 5 % full donor allele (blocking-only conversion) where the tool reports > 10 %.
5. **Allele recovery.** For each truth allele (indel or two-cut deletion, summed over alleles of equal net size) at >= 10 %: recovered if the tool lists that net size at between half and twice the true fraction. Reported overall and for large deletions (40, 80, 150 bp).
6. **KO-score error.** MAE of the frameshift-or-large share (alleles with size not a multiple of 3, or >= 21 bp) for tools that report one.
7. **Coverage and speed.** Share of traces analysed; median seconds per trace on this host.
8. **Per-family table.** MAE for each of 15 scenario families, so a gain in one family cannot hide a loss in another.

## Pre-registered success criteria for "beats the comparators"
- Primary: on the test loci, Assured QC's edited-fraction MAE (analysed traces) is lower than that of every comparator, and the 95 % bootstrap interval of the paired per-trace difference excludes zero.
- No-regression rule: in no scenario family is Assured QC's MAE more than 2 points above the best comparator's.
- Safety: false-positive share on wild-type replicates <= 2 %; false-intended rate on blocking-only conversions = 0; coverage >= 99 %.
- Targets (reported whether or not met): MAE <= 2 points for truth >= 10 % and <= 1 point for truth <= 5 %; detection limit <= 2 %; allele recovery >= 95 % for alleles >= 10 %; 95 % interval coverage 90-98 % (step 6).
- Anything not met is reported as not met.

## Known limits of the evidence
D1 truth is exact but the traces are generated (limits listed in `bench/simlib.py`); D3 is a protocol, not data; no comparator was tuned by us, and ICE is run with its published defaults on the same inputs; web-only tools are absent from every table.
