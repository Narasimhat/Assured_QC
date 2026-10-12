# Marker-level peak-height response: what was tested and why no correction was added

Question: do donor markers read with a marker-specific peak-height response (an alternate base that never reaches 100% of the signal, or a 50:50 mixture that reads 41% at one marker and 63% at another),
and would a per-marker correction improve the intended-edit estimate?

## What the simulator can and cannot say

The held-out simulator (`bench/simulate.py`) draws peak heights from real control traces and applies dye scaling of up to 10%; it contains no marker-specific response of the size seen in the real
clones below. A fresh simulation seed therefore cannot test the hypothesis: a model tuned on it would be fitted to the noise the simulator imposes. The fresh-seed held-out run was not done for this item.

## What the real clones show

`bench/marker_dump.mjs` writes the background-corrected alternate-base share at every marker for each SNP clone. Four of the six SNP sets have two or more markers (45 clones in all):

- In one set (10 clones), 7 clones read 0.99 at the first marker and 0.78 at the second, identical to two digits. That is a marker-specific ceiling, not clone-to-clone variation.
- Estimating each marker's ceiling from the other clones of the batch (median of the shares at or above 0.65, leaving the clone out) and rescaling shrinks the mean between-marker spread in that set from 0.234 to 0.067
  and the distance to the nearest of 0, 50 and 100% (the clonality prior) from 0.098 to 0.013.
- In the other three sets (13, 11 and 11 clones) the ceilings are 0.97 or higher and the correction changes nothing (spread 0.265 to 0.265, 0.861 to 0.862, 0.149 to 0.143).
- In the one set where the ceiling matters, the allele-level fit already reports an intended edit of 99.7% or more for those 7 clones; the engine's expected composition at a marker comes from the control's own peak pattern,
  which absorbs a response of this kind. The marker-level readout (0.78) is not the quantity the intended-edit estimate uses.

## Decision

No per-marker response model was added to the estimator. The effect appears at one locus, the engine's allele fit already handles it there, and no truth exists to show an improvement. The real-data guard is unchanged
by the work in this release: knockout clones differ from the chromatogram readout by 2.88 points on average (68 clones) and 3 of 84 SNP clones read an intended edit more than 10 points above the donor-base bound.
The weaknesses listed in `docs/results.md` (wild-type false positives, a detection limit of about 5%) come from allele-level noise absorption and are not addressed by marker ceilings.

What would settle it: the wet-lab known-mixture plates (`docs/known-mixture-protocol.md`) give marker readouts with known fractions, from which a marker's response can be measured with truth and a correction tested on a held-out plate.
