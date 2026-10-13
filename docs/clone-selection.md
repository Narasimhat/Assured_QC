# Clone selection

`src/selectionFigure.js` draws one figure from the design to the clones to take forward, and the report gets a section "Clone selection" with the same figure and the suggested next steps.
Browser: a section under the summary graph (SVG and PNG, 300 dpi). Command line: `clone_selection.svg`. Vector SVG with real text, no scripts and no external resources.

## What is in the figure

- **a. Clones to carry forward first, aligned with the wild-type control.** The chromatograms of the clones named by the recommendation, under the guides (grey bar, PAM darker), cuts and designed changes of the design (`docs/aligned-chromatograms.md`).
- **b. All clones.** One row per clone, grouped by what to do with it (carry forward first, reserve, backup, sequence again, not recommended, discard; four or more discards are one row): the share of the donor base at each designed change from the peak heights, the genotype read from the sequence, the call (A accept, H hold, S sequence again, R reject) and the support (the lower bound of the intended-edit interval for an edit, the knockout support for a knockout). A clone with a quarter or more of its signal in indel alleles shows n/a at the designed changes, because its bases there are shifted.

## How many clones to carry forward

A knockout takes one clone per genotype, so three by default. For a wanted edit (SNP, small tag, knock-in) clones of one genotype are replicates and all count, so the default is six. `--panel N` (command line) or the box under the summary graph (browser) sets it; the ranking does not change, only where the line between "carry forward first" and "reserve" falls. Clones whose sequence does not support the wanted edit are never in the panel.

## Next steps

For the clones carried forward, the report lists suggested steps, which a laboratory can adapt: for an edit, validation in both directions with a second PCR outside the first primer pair (allelic dropout, a large deletion on one allele), phase of silent changes that differ between alleles, off-target sites; SNP-array karyotype; banking; quality control (mycoplasma, STR identity, pluripotency markers, whole-genome sequencing of the final clones). For a knockout, protein level, a PCR across the cuts, karyotype, banking and quality control.

## Design reports with several donors

An exported design report lists the donors to order. For a point edit, the ssODNs recommended for each guide are used together: the markers of all donors are merged, the change every donor carries is the intended edit, a change some donors carry is a blocking change, and each donor lists the markers it carries (`carriesMarkers`), so a clone made from either donor is read as the intended edit and a clone made from both is read as a compound heterozygote. Before 0.5.0 only the first recommended donor was read from a report, so the changes of the other donors were not known and clones made from them were called partial conversions. Donors with an insert or a replacement are not combined (the first donor is used and a note says so); `--donor SEQ1,SEQ2` gives several on the command line.
