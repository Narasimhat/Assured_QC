# Changes in 0.3.0 (preview) and the evidence for them

0.3.0 adds two small, separately switchable checks to 0.2.0. The allele library, the fit, the intervals and the decision rules are unchanged, and on 350 simulated traces the output is identical to 0.2.0. The checks were chosen by running Assured QC next to two other engines (sangerdecon, a Python ICE-like tool; Assured_TraceEdit) on the lab's own pairs and on a fresh simulated set, and porting only what the real pairs showed Assured QC missing.

## 1. Offset re-anchoring at the cut (`src/offset.js`, option `refineOffset`)

**Problem.** The edited read is laid onto the control with one offset, measured on the high-quality stretch upstream of the cut. If one of the two reads has a base call more or less between that stretch and the cut (a merged or split peak), every position after the step is compared one base out. On one real control read that carried two extra calls 20-60 bases before the cut, all 14 clones of a two-guide knockout project were reported one base too large (31, 47 and 49 instead of 30, 46 and 48 bases) and the one unedited clone was reported as a one-base deletion at 100 %. Aligning the called bases of each clone read to the control gave clean 30-, 46- and 48-base gaps with no gap in the read, and the same two steps in the offset in all 14 clones, including the unedited one, so the steps belong to the control read, not to the alleles. ICE 1.2.0 reports the same one-base-larger size on a clone of that project.

**Check.** In the 20 bases that end 2 bases before the first cut (where an edit cannot yet have started) the called bases of the two reads must agree under the offset in use. Only if at most half of them agree, and a shift of at most 3 bases makes at least 85 % agree over at least 15 positions, is the offset changed; the old offset still applies to positions before the point where the new one starts to agree, so the upstream stretch is compared as before. A pair whose offset already fits is never touched. A warning states the shift.

**Result on the real pairs.** 18 of 195 pairs changed. In the project above 11 of the 13 deletion clones now have the correct size (R2 0.90 to 0.998) and the unedited clone is now wild type. In a second project (7 clones, low-quality reads) four clones and the bulk pool changed, and all 7 clones now agree with sangerdecon's allele sizes (3 did before). Across the 174 clones analysed by both Assured QC and sangerdecon the allele sets agree exactly in 148 (132 with 0.2.0) and within one base in 151 (149).

**Limits.** Two clones of the first project, whose deletions start 10-12 bases upstream of the cut, stay one base off because the deletion itself lies inside the 20-base check window. A drift inside the last ~22 bases before the cut cannot be told from an edit and is not corrected (tested). A real 1-3 base deletion that starts more than ~15 bases upstream of the cut and does not cross it would be absorbed into the offset; this is rare for Cas9 and was not seen. The real-pair evidence is one control read.

## 2. Unplanned base changes (`src/unplanned.js`, option `unplanned`)

**Problem.** The library holds what the design says can happen. A clone that also carries a base change the design did not ask for is fitted as well as the library allows and shows only as a lower R2. In a SNP project four clones carried a heterozygous T>C ten bases upstream of the cut (raw peaks: a second base at 48-52 % there, 0 % in the other nine clones); 0.2.0 reported them as R2 0.88-0.89 against 0.99+ for most of the other clones.

**Check.** A position is reported when a base other than the control's carries at least 25 % of the edited signal, the fit predicts less than 10 % of that base there, the control is clean at that position, and it is not a planned marker. Positions upstream of the cut are always scanned; positions after the cut only when every allele in the fit leaves the length unchanged. More than two hits is a misfit and nothing is reported. Positions within 3 bases of a re-anchoring step are skipped. The result carries `unplanned` (position relative to the cut, change, share), a warning, a CSV column `unplanned_changes`, and a clone that would otherwise be accepted is held with the reason (`holdOnUnplanned: false` turns that off). Numbers, intervals and genotype are not changed.

**Result.** The four clones are named (-10, T>C, 49-52 % in the fit). Across the 195 real pairs 10 pairs are flagged: the four, three more in which raw peaks also show a second base at a similar share at that position, and three whose position could not be mapped to a base call (unverified). The cause of the three confirmed extra hits is not known (a shared variant, an allele outside the library, or a trace artefact). On the 348 analysed simulated traces of the second fresh set, which contain no unplanned change by construction, one trace was flagged (a noisy scenario): 0.3 %.

**Limits.** Isolated positions only. It names a position, not a cause.

## Evidence set up before the results were seen

`bench/PREREGISTRATION_0.3.0.md` was written before the two fresh seeds were generated. Seed 20261010 (the 0.2.0 README seed) was used only to check that the pipeline reproduces; it was not used for any decision.

**Fresh seed 20261201 (350 traces, 12 loci, forward read only), edited-share error in points (95 % CI), frozen `bench/score.py`:**

| Tool | MAE | WT false positives (>= 2 points) |
|---|---|---|
| Assured QC 0.2.0 | 1.94 (1.56-2.39) | 6 % |
| sangerdecon 0.1.0 | 3.15 (2.47-3.91) | 11 % (3 % with its own detection limits) |
| Assured_TraceEdit | 3.56 (2.75-4.46) | 14 % |
| ICE 1.2.0 | 6.21 (4.73-7.84) | 3 % |

Paired differences over common traces (positive = Assured QC better): sangerdecon +1.02 (0.42 to 1.72), TraceEdit +1.62 (0.98 to 2.34), ICE +4.17 (2.70 to 5.82). Assured QC has the lowest error in 6 of 14 trace families and is within 1.0 point of the lowest in 13; no family met the pre-registered rule (another tool better by more than 1.0 point, paired CI excluding 0, at least 15 traces), so nothing was ported on accuracy grounds. The loci are the lab's own: Assured QC and sangerdecon were developed on controls from the same data, so "fresh" means new mixtures, not new loci.

**Fresh seed 20261215 (350 traces), 0.3.0 against 0.2.0 (pre-registered criteria):**

| Criterion | 0.2.0 | 0.3.0 | Met |
|---|---|---|---|
| Edited-share MAE, analysed (95 % CI), not worse by more than 0.3 | 1.84 (1.49-2.28) | 1.84 (1.49-2.28) | yes |
| WT false-positive rate not higher | 3 % | 3 % | yes |
| Coverage not lower | 99 % | 99 % | yes |
| New features fire on synthetic positives and stay silent on negatives | | `test/offset-and-unplanned.test.js` | yes |
| All earlier tests pass | 107 / 107 | 115 / 115 (107 + 8 new) | yes |

All 350 outputs are identical between versions, and the re-anchoring never triggered on a simulated trace (0 of 348): the simulator renders each edit from the control's own peaks, so it contains no base-caller drift. The simulated sets therefore show that nothing regressed, not that the checks help; the help is shown by the real pairs above.

## Not done

* The browser UI was not exercised by hand. The production build (`vite build`) and the engine tests, including the CLI test that cannot run in the author's sandbox, pass in CI on this branch.
* The README criteria 0.2.0 did not meet (WT false positives, 5 % detection limit, graded-HDR error, interval coverage) are unchanged.
* ICE 1.2.0 (non-commercial research licence) was run locally as a comparator and is not included here.
* Merging remains held for the owner's decision; this branch is stacked on `feat/assured-qc-0.2.0`.
