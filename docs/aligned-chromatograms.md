# Clones aligned with the wild-type control

The decision and the allele fit say what each clone is made of. These figures show it on the traces: every clone laid onto the wild-type control read, base by base, with the guides, PAMs, cuts and designed changes marked.
`src/alignedView.js` builds the view of one sample (plain data, kept on the result as `view`), and `src/alignedFigure.js` draws it. Vector SVG with real text, 1 unit = 1 pixel at 96 dpi, no scripts and no external resources.

## What is drawn

| Figure | Where | What |
| --- | --- | --- |
| Chromatogram aligned with the wild-type control | one per sample in the HTML report, in the browser sample view, `chromatograms/<sample>.svg` with `--clone-figures` | The control chromatogram above the clone's, on the same stretch. Guides as grey bars (PAM darker, arrow towards the PAM), cuts as dashed red lines, designed changes as shaded columns (orange: blocks a PAM, purple: the edit). The base calls of both reads above their curves; a base that differs from the control is boxed in red, and a second base that is at least 20% of the signal is printed small above the call. |
| All clones aligned, base calls | report section "All clones aligned", browser, `clones_alignment.svg` | One row per clone, best clone first (the order of the run summary), on the control sequence: a dot where the clone agrees, the coloured base where it differs (half-and-half or dashed where two bases share the position), a dash where it lacks the base, a purple triangle where it has extra bases. Rows for hundreds of clones stay readable. |
| All clones aligned, chromatograms | the same | The control and every clone's chromatogram stacked on the same stretch, up to 48 clones (`clones_aligned.svg`). |

Clones are grouped by design: samples with the same guides, designed changes and stretch share a figure; a run with several designs gets one figure per group (`clones_alignment_2.svg`, ...).

## How a clone is laid onto the control

The stretch runs from 18 bases before the first guide, PAM, cut or designed change to 28 bases after the last one (40 to 220 bases). The clone read is mapped to the control with the offset the analysis used (`prepare.js`, re-anchored at the cut when the base callers disagreed upstream).
Only when a mapping with gaps explains the clone's base calls clearly better is it drawn: an affine-gap alignment of the clone's calls to the control's (match 2, mismatch -2, gap open -8, extend -1) is used when it has at most three gap or insertion runs, at least 90% identity over the stretch, and at least 15 points more identity than the offset mapping. A homozygous deletion or insertion is then drawn as a grey block of the missing bases (or a purple triangle for extra bases) and the bases after it line up with the control again. A clone with a wild-type or different allele, or a noisy read, keeps the offset mapping, so the second allele shows as mixed signal after the cut instead of being forced into a tidy picture.

Curves are the raw trace, resampled so that every control base is one unit wide (6 points per base, peak locations as the base centres) and scaled to the tallest peak of the stretch, control and clone separately.

## Limits

- Only the first read of a clone is drawn; a clone read from both primers is fitted jointly but shown on its first read. A read from the reverse primer is drawn in its own direction (the design is turned around, as in the analysis).
- The figure shows what the base caller and the trace say. The allele shares and intervals come from the fit, not from the picture; a second base is flagged from the composition of the read at that position and is not a call of heterozygosity.
- Per-sample figures are drawn for every sample of a run of up to 48, and for samples that were not rejected above that. The stacked chromatograms stop at 48 clones; the base-call alignment has no limit.
- `results.json` leaves the view out (as it leaves the charts); the figures need the reads, so they are made by a run and not by `recommend`.
