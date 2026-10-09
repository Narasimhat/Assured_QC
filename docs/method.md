# Method

## Input and orientation
Each trace is parsed from the ABIF file (called bases, peak locations, four channels, quality values). The design reference is the amplicon sequence in the design file, or the control read itself when guides and donor are typed in. If the read matches the reverse complement better, the design is turned around, not the trace.

## Alignment and window
The control and the edited read are aligned to the reference with a local Smith-Waterman alignment. The offset between the two traces is taken from a local alignment of the edited calls to the control's longest high-quality stretch before the cut (the diagonal most bases near its downstream end lie on). If the edited read's first calls are poor and only part of the stretch aligns, the offset is accepted when at least 30 bases (40% of the stretch) align at 95% identity or better, and a warning says so. Positions where the control already differs from the reference (parental SNPs) are patched in the reference. The inference window runs from 25 bp before the first cut to 150 bp after the last cut.

## Peak model
At each position the four channel heights around the called peak are normalised to a composition (fractions of A, C, G, T). The control composition at the same context gives the background pattern of the sequencing reaction (shoulder peaks, dye-specific heights).

## Allele library
Candidate alleles are built from the reference: wild type; every indel of 1 to 30 bp at each cut with every placement; +1 and +2 insertions; the deletion between cuts; the intended edit (every donor difference from the reference); partial conversions (contiguous runs of donor changes and all-but-one-marker sets). Alleles with identical sequence are merged and reported with the alternatives they cannot be told from.

## Deconvolution
The predicted composition of a mixture is the weighted sum of the allele compositions; weights are found by non-negative least squares (Lawson-Hanson), with a flat background column. Where an allele's local context equals the wild-type context, the control composition is used for it, which removes sequencing artefacts that would otherwise be read as editing. Only the generic indel family is shrunk towards zero by the noise level; wild type, the intended edit and partial conversions are not shrunk. Small penalties order otherwise tied explanations (wild type and edit first, partial conversions and skipped features next, generic indels last). R-squared is reported relative to the wild-type-only fit, so it carries little meaning when nothing is edited.

## Summaries
Editing % is 100 minus wild type. Signal that no allele in the library explains (the background term of the fit) above a 10% baseline is reported as "unexplained" and counted as not wild type; the allele percentages are then shares of the total signal and sum to 100 minus the unexplained part. Without this, a complex pool of many rare alleles would be renormalised onto the few alleles that were fitted and read as wild type. Knockout score is the fraction of alleles with an indel that is not a multiple of 3 or is 21 bp or longer. Intended edit and blocking-only conversion are separate rows. Clone calls: homozygous if one allele is at least 85%; two alleles each at least 25% and together at least 85% give the heterozygous or compound categories; anything else is "mixed".

## Junction reads
A read is aligned to the expected complete knock-in allele (genome flank, arm, insert, arm, genome flank). The report lists differences by region, whether designed blocking changes are present, which junctions are spanned, and how far the read extends into genome flank. If the wild-type sequence explains the read better than the knock-in allele, the read is refused as evidence for the knock-in.

## Out-out band sizes
Observed sizes are matched (within 5% or 20 bp) to the wild-type, knock-in and deletion amplicon sizes in the design.
