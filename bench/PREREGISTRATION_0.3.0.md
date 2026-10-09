# Pre-registration — three-tool comparison and port decision (0.3.0)

Written 2026-10-09 BEFORE generating or scoring seeds 20261201 and 20261215. Nothing below may be changed after those results are seen; deviations get listed in the report under "Deviations".

## Data
* Simulator: `bench/simulate.py` from Assured_QC PR #1 (frozen), run on one real control per analysable project from the lab's own data (twelve loci; names withheld, as in the 0.2.0 public copy). Projects whose control cannot host the scenario geometry are skipped by the simulator and reported.
* Seed 20261010 (the seed published in the 0.2.0 README) is used ONLY to check that the pipeline reproduces; it is not evidence for any decision. Seed 20261201 = decision set (Phase 2). Seed 20261215 = no-regression set (Phase 3). No code is tuned on either.
* Limits stated up front: the loci are the lab's own, and Assured QC and sangerdecon were each developed on controls from these same data. "Held out" therefore means held-out mixtures (new seed, new noise, new fractions), NOT held-out loci. Absolute errors are optimistic for every tool in the same way; the between-tool differences are the finding.
* Simulated traces are rendered from real control peaks; this favours engines that use the control's own peak shapes (all three do). PCR length bias, independent reverse runs and novel-junction context effects are not modelled.

## Tools and settings
Assured QC (PR #1 branch, default options, workflow chosen by the harness), Assured_TraceEdit engine (deletion ≤ 40, insertion ≤ 2, as in the harness adapter), sangerdecon 0.1.0 (default options; two variants: raw point estimate and with its own detection limits applied). Peaksplit and ICE are scored only if they install here; ICE is a comparator only (non-commercial licence). Forward read only for all tools (what the frozen adapters pass).

## Metrics (frozen `bench/score.py`, `--split all`; definitions in docs/scorecard.md)
Edited-share MAE on analysed traces (bootstrap 95 % CI), MAE with failures imputed as 0 edited, intended-edit MAE, KO-score MAE, WT false-positive rate (WT traces with estimate ≥ 2.0 points), per-family MAE, large-deletion recovery, sensitivity by level (detection limit). Paired differences vs Assured QC use the scorer's paired bootstrap over traces both tools analysed.

## Decision rules for ports (Phase 3)
1. A capability is ported from another tool when EITHER (a) on the real lab pairs it surfaces a finding that Assured QC's output does not contain and that is independently supported (e.g. a control-vs-sample base change confirmed in both reads, or a genotype contradicted by a ploidy impossibility), OR (b) on seed 20261201 another tool beats Assured QC in a trace family by more than 1.0 point with a paired 95 % CI that excludes 0 and n ≥ 15 traces in that family.
2. Nothing is ported because it is "also nice to have". A family win by a tool that fails on more traces than Assured QC is judged on the common set and the coverage difference is reported.
3. If the evidence shows that a different scope is needed (e.g. a different engine would have to replace Assured QC's), stop and ask the owner instead of building.

## No-regression criteria for 0.3.0 vs 0.2.0 (seed 20261215)
* Edited-share MAE (all, analysed) not worse than 0.2.0 by more than 0.3 points.
* WT false-positive rate not higher than 0.2.0.
* Coverage (fraction analysed) not lower.
* Every new feature fires on its synthetic positive fixture and stays silent on its synthetic negative fixture (unit tests).
* All 107 existing tests still pass.
Anything unmet is reported as unmet.
