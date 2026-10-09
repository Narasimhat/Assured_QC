# Validation

The held-out simulated scorecard against ICE, peaksplit and Assured_TraceEdit, with the pre-registered criteria and what was not met, is in `docs/results.md`; the protocol and success criteria were fixed before it was run (`docs/scorecard.md`, `bench/FROZEN.json`). The sections below predate it.

All numbers are percentage points. "Truth" is the known mixture used to generate the synthetic trace. ICE is Synthego ICE 1.2.0 run locally.

## Synthetic traces (17 pairs, `fixtures/sim`, `docs/validation-synthetic.csv`)

Mean absolute error against the known mixture (largest single error in brackets):

| Quantity | This engine | ICE 1.2.0 |
| --- | --- | --- |
| Editing % (100 - wild type) | 1.5 (4.5) | 1.6 (7.0) |
| Intended edit (HDR) | 1.2 (5.9) | 4.0 (43) |
| Knockout score | 1.1 (2.7) | 0.5 (3.0) |

The ICE HDR error is dominated by one clone with 50% intended edit and 50% blocking-changes-only alleles (ICE reports 93%) and a pool with 35% intended edit plus 10% blocking-only (ICE reports 50%). ICE counts every read that carries the donor-derived changes as HDR; this engine separates them. The synthetic traces were used while this engine was built and tuned, so the engine's figures are optimistic compared with unseen data.

## ICE example traces (from the ICE repository; the traces are not redistributed)

| Example | This engine | ICE |
| --- | --- | --- |
| Good example (single guide), editing % | 73 | 77 |
| +14 knock-in, HDR % | 35 | 33 |
| +14 knock-in, editing % | 92 | 90 |
| +14 knock-in, knockout score | 39 | 39 |
| 2-bp substitution, HDR % | 34 | 44 |
| Three-guide multiplex, editing % | 74 | 63 (ICE-d 61) |

The 2-bp substitution differs by 10 points; reading the two substituted positions directly from the peak heights gives 36 to 38%, closer to this engine. The three-guide multiplex example does not agree and is a known weakness; use it with care.

## Laboratory traces: 172 control/edited pairs against ICE and against the chromatograms (clone-level table not published)

Data: Sanger traces from this laboratory, 12 edit sets in hiPSC: knockouts with two guides (L01, L02, L03, L04, L05 exon 4, L06; 68 clones) and SNP corrections with one guide and an ssODN (L07, L08, L09 exon 20, L10, L05 exon 25, L11; 92 clones), plus one bulk pool per set (12). Each edited trace was paired with the wild-type control of its set. Both tools received the same control, edited trace, guides and ssODN. ICE 1.2.0 was run locally with default settings (a Biopython compatibility patch was applied at run time). The traces are not redistributed.

**There is no sequencing-confirmed truth for these samples.** Accuracy was therefore judged against three checks that use the chromatograms directly and none of the tools' models (all three share the trace noise with the tools; they are checks, not a gold standard):

1. *Knockouts, edited share.* In the 100 bases after the first cut, the share of the signal that sits on the wild-type base (control, aligned in phase) is converted to a wild-type fraction f = (share / pre-cut purity - 0.25) / 0.75, the 0.25 being the chance that another allele shows the same base. Edited share = 100 (1 - f). It is blind to a 1 bp indel flush with the cut.
2. *SNP corrections, donor bases.* At every position where the ssODN differs from the control read, the share of the donor base against the control base in the edited trace. The smallest of these shares is an upper bound for the fraction of alleles that carry the complete donor change; an "intended edit" above it (plus 10 points) is not possible. Alleles that carry the donor bases and an indel are not intended edits, so a value well below the bound is allowed.
3. *Clonality.* A clone has 0, 50 or 100% edited alleles. The distance of each estimate from the nearest of these three values measures noise (mosaic clones would count against both tools).

| Check | Assured QC | ICE 1.2.0 |
| --- | --- | --- |
| Knockout clones (n = 67 analysed by both): mean absolute difference from trace readout, edited share | 2.9 (r 0.97, mean bias +1.9, 1 clone off by more than 15) | 14.3 (r 0.72, mean bias -11.4, 24 clones off by more than 15) |
| SNP clones with a readable donor position (n = 84): intended edit more than 10 points above what the donor-base readout allows | 3 | 17 |
| SNP clones with a clean trace (n = 66): intended edit more than 10 points below the readout | 2 | 11 |
| All clones analysed by both (n = 159): estimate within 10 points of 0, 50 or 100% | 95% | 67% |
| Bulk pools, knockout (n = 6): mean absolute difference from trace readout | 3.9 | 10.1 |
| Bulk pools, SNP (n = 6): mean absolute difference of intended edit from the donor-base bound | 1.6 | 19.5 |
| Median time per pair (one core) | 0.2 s | 10.2 s |

The 10-point tolerance covers the error of the readout itself and is the same for both tools. At zero tolerance Assured QC exceeds the bound in 20 clones (14 of them by 1 point or less) and ICE in 17 (all by more than 10 points).

Direct comparison of the two tools (clones): knockout edited share r 0.77, Assured QC higher by 13 points on average; knockout score r 0.91, mean difference 0.2 points; SNP edited share r 0.84; SNP intended edit / ICE HDR r 0.77, mean absolute difference 14.5 (29 of 92 clones differ by more than 15 points).

What the differences are:

- **Blocking-only conversions counted as HDR by ICE.** In the L09 exon 20 set (13 clones and the pool) the donor base at the intended position is absent from every allele (readout 0) while the other two donor changes are present on all alleles in 10 clones. Assured QC calls these ten "homozygous_partial" with intended edit 0.6-1.0% in nine and 7.7% in one; ICE reports HDR above the readout bound in 12 of the 13 clones (68% in nine, 37-70% overall) and 62% for the pool. This is the case the synthetic traces predicted. In three L08 clones ICE HDR exceeds the readout bound by 14-24 points where Assured QC separates 43-57% blocking-only alleles.
- **Complete conversions under-reported by ICE.** Nine L11 clones have the donor base on every allele (readout 100%) and Assured QC reports 99.5-100% intended edit; ICE reports 84% for all nine. The cause was not determined.
- **Two-guide knockouts.** For clones whose downstream signal is a single sequence that differs from the wild type (38 clones, ICE result available), Assured QC reports 99% edited on average and ICE 85%. Eleven of the 14 L04 clones are called as the same homozygous 31 bp deletion between the two cuts (Assured QC 100%, ICE 80% in all eleven); one has a 1 bp deletion and two traces are unusable.
- **ICE gave no result for one clone** (quality too low); Assured QC analysed it with warnings.

Genotype calls (160 clones): homozygous edit 43, homozygous partial 10, homozygous indel 26, edit plus indel 17, compound heterozygous 13, heterozygous indel 5, heterozygous edit 2, edit plus partial 3, wild type 9, mixed 29, unclear 3. All 43 "homozygous edit" calls have the donor base on every allele in the readout, all 10 "homozygous partial" calls have none at the SNP position, and all 4 "homozygous wild type" calls with a readable donor position have none.

To repeat the comparison on your own pairs: `scripts/run_validation.mjs` (this engine), `scripts/run_ice_jobs.py` (a harness that calls a local ICE installation; ICE is not included), `scripts/readout_validation.py` (the chromatogram readouts). Each takes a JSON list of control/edited pairs with guides and optional ssODN.

### Failures found on these data and fixed

These two failures were found by running the engine on this data set, so the numbers above are no longer blind for them.

1. A read whose first calls are poor was refused ("does not align", 1 of 172 pairs): the upstream stretch of the control was required to align over 60% of its length. It is now accepted when at least 30 bases (at least 40% of the stretch) align at 95% identity or better, with a warning. Regression test: `test/real-data-regression.test.js`.
2. A complex pool was reported as 4% edited (trace readout 87%, ICE 76%). The fit assigned 94% of the signal to background and the 6% left was renormalised to 100%, so the one wild-type column dominated. Signal that no allele explains above a 10% baseline is now reported as `unexplained` and counted as not wild type (84% edited, 83.5% unexplained for this pool), with a warning from 15% upwards. Regression tests in the same file. For knockout clones the correction lowered the mean difference from the readout from 3.4 to 2.9 points and left the nine clean wild-type clones below 2%.

### Known weaknesses on these data

- Three L05 exon 25 clones with poor traces get an intended edit of 16-34% where the donor base is absent (two have pre-cut discordance 0.58-0.59, the third a fit R-squared of 0.33); the engine warns that the fractions are unreliable for all three. Poor traces are the main source of error for both tools.
- Warnings are frequent: 94 of 172 pairs carry at least one reliability warning and 78 the low-quality warning. For knockout clones the error does not depend on these warnings (mean difference 2.8 with, 3.7 without the low-quality warning); for SNP intended edits it does (10.4 with, 6.3 without).
- Compound heterozygous and mixed knockout clones are called but not resolved allele by allele: 29 clones are "mixed".
- The trace readouts above cannot see a 1 bp indel flush with the cut, and the donor-base bound says nothing about indel alleles.

## Failure cases covered by tests (`test/safety.test.js`)

- The control analysed as its own edited sample: 0% editing, homozygous wild type.
- Control and edited traces swapped: refused (control does not match the design reference).
- Reads that stop before or just after the cut: refused with the reason.
- Noise added to the edited trace: at a noise standard deviation of 50 (peak height about 2,100) the intended edit is still within 6 points of the true 50%; at 250 the pre-cut discordance exceeds 0.25 and a warning is raised. At 400 the unwarned result would have been 19% intended edit instead of 50%, which is why the warning exists.
- A damaged .ab1 in a batch is reported by file name and the other samples are analysed.

## Not yet validated
Clones whose genotype is confirmed by sequencing (amplicon NGS or cloned Sanger). The laboratory data above were checked against ICE and against the chromatograms, not against a confirmed genotype, so treat the output as a second opinion next to ICE until a subset of clones has been confirmed.
