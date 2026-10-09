# Design hand-off and expected alleles

The QC app reads one JSON object per design (`assured-qc-design/1`), exported by the design app
(`src/qcExport.js`, `buildQcDesignSpec(result)`). Everything the analysis needs is in it, so the QC app does
not depend on the design engine.

> The design app's exported HTML report is also a hand-off: Quick Start and the classic form read its guides and recommended ssODN donor (`src/designReport.js`). A report has no reference window, so the junction check for large knock-ins still needs the JSON spec described below.

## Coordinates

All positions are 0-based indices into `reference.sequence`, a window of the uploaded reference (the + strand)
that starts at `reference.offset`. A cut index `c` means the cut lies between bases `c-1` and `c`.

## Fields

| Field | Content |
|---|---|
| `design` | `type` (ko, pm, it, ct, nt), `editKind` (deletion, snp, internal_tag, tag_small, tag_large, reporter), gene, one-line description of the intended edit |
| `reference` | window sequence, its offset in the uploaded reference |
| `guides[]` | name, spacer, PAM, strand, `cut`, protospacer start/end on the window |
| `donors[]` | name, format (ssodn, block, aav), + strand sequence, `refStart`/`refEnd` (the reference span it replaces), `insertBp`, `replacedBp` (native codon removed by the cassette), arm lengths, `insertStart`, matched guide, ssODN order strands |
| `markers[]` | every base the donor changes outside the insert: `pos`, `ref`, `alt`, `role` (`intended` = the SNP; `blocking` = silent or PAM change), the guide it blocks, label |
| `primers[]` | name, sequence, placement on the window |
| `amplicons[]` | primer pair, wild-type size, edited size per donor (and the deletion size for two-guide designs) |

`validateDesignSpec` (`src/designSpec.js`) refuses a hand-off that disagrees with itself: spacer not at its stated
position, cut not 3 bp from the PAM, marker base not the reference base, a donor difference that no marker
declares, primers not on the window, amplicon sizes that do not follow from the primers.

## Expected alleles

`buildExpectedAlleles(spec)` returns every sequence a pool or clone can plausibly carry, each with the full
window length:

* **wt**
* **indel**: deletions of 1-30 bp that remove the cut (every placement from "starts at the cut" to "ends at the cut"),
  and +1/+2 insertions with every possible inserted base, at each guide.
* **deletion_between_cuts**: loss of the stretch between two cuts, with end positions varied by up to 3 bp.
* **edit / edit_partial**: the donor as a run of features in reference order (each marker, and the insert if there is one).
  Homology-directed repair copies a contiguous run, so every contiguous run is an allele; the full run is the intended edit.
  A single marker can also be lost from an otherwise complete run (that mismatch repaired back), so every "all but one
  marker" set is included. The insert is never lost this way. Partial alleles are what let the app say "blocking changes
  copied but the SNP not" instead of counting it as a knock-in.

Alleles with identical sequences are merged; the more specific description wins (edit > partial > deletion between cuts > indel)
and the other is kept as an alias.

## Why partial alleles matter (measured)

Synthego ICE 1.2.0 on synthetic APOE R176C traces (fixtures/oracle_vs_truth.csv): a clone that is 50% intended edit and
50% blocking-changes-only (SNP not copied) is reported as 93% HDR, because the donor-derived changes dominate the fit. A pool
with 35% intended and 10% blocking-only edits is reported as 50% HDR.

These traces are synthetic (linear mixture model, own generator); behaviour on real traces is to be confirmed on your own
control/edited pairs in the validation step.
