# Pre-cut discordance warning

The warning "The control and edited traces already differ before the cut" compares the control and the edited read over the stretch upstream of the first cut,
where both carry the same sequence. It is meant to catch a noisy read or a control from a different template.

## What changed

- A knockout or deletion design is read from the cut and what lies after it, so for these designs the warning now looks at the 40 bases next to the first cut
  (`NEAR_CUT_BASES`, `quality.discordanceNearCut`). A mismatch early in the read, where the stretch next to the cut still matches (typically a base-calling shift
  in the noisy first part of a read), no longer warns.
- A design that reads donor markers (SNP, small tag, base editing) keeps the whole alignment window, because the bases it reads can lie anywhere in it. It also
  warns when the stretch next to the cut differs.
- Estimates do not change: the maximum difference in the edited share over 160 laboratory clones is 0.0 points. Only the warning and the line in the report change
  (the report now prints the figure for the stretch next to the cut beside the whole-window figure).

## Evidence (160 laboratory clones: 68 knockout or deletion, 92 SNP with ssODN)

There is no sequencing-confirmed truth for these clones. Error is measured against the model-free chromatogram readouts of `docs/validation.md`: the wild-type base's
share of signal after the cut for knockouts (error = engine edited share minus readout) and the donor-base share at every marker for SNP clones (an overcall is an intended
edit more than 10 points above the readout's upper bound).

| | whole-window rule (before) | rule now |
|---|---|---|
| warned, all 160 clones | 46 | 29 |
| warned, 68 knockout or deletion clones | 39 (57%) | 20 (29%) |
| warned, 92 SNP clones | 7 | 9 |
| donor-base overcalls (3 SNP clones) that were warned | 3 of 3 | 3 of 3 |
| knockout clones more than 5 points from the readout (11) that were warned | 4 of 11 | 1 of 11 |

For knockouts the whole-window figure did not predict error (AUC 0.49 for error above 5 points; mean error 3.1 in the clones it flagged against 2.6 in the others), and the
figure next to the cut does not either (AUC 0.40). The warning therefore has no demonstrated value for knockout accuracy; the change removes alarms that
did not carry information, it does not add a predictor. For SNP clones the whole-window rule flagged all three overcalls, which is why it is kept for designs that read
markers. Caveats: three events, all from one locus, and the readouts are not sequencing truth. Treat the SNP rule as retained, not as shown to work.

Re-run: build a jobs file as in `scripts/run_validation.mjs` and run `node bench/discordance_eval.mjs jobs.json out.json`; the output lists both discordance figures and the
warning per clone.
