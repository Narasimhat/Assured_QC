# Known-mixture protocol (wet-lab ground truth for Assured QC)

Purpose: obtain Sanger traces whose allele fractions are known from the bench rather than from a simulation or from another tool. These traces are the physical anchor of the benchmark in `docs/scorecard.md`. Nothing in this protocol has been run yet; the numbers in `known-mixture-*.csv` are the design.

## 1. Design in one paragraph
Amplicons from the parental line (wild type) and from genotype-confirmed homozygous clones are purified, quantified, and mixed **after PCR** at defined molar fractions (0, 2, 5, 10, 25, 50, 75, 100 %), so PCR length bias cannot distort the truth. Every mixture is sequenced from both ends in the same run, two independent preparations per level, with run replicates for eight mixtures. Analysts receive only trace files and a coded plate layout; the key stays sealed until the analysis settings are frozen.

## 2. Series (74 mixtures, 166 reactions on 2 plates)
| Series | Locus / question | Mutant source | Levels (mutant %) | Preparations |
|---|---|---|---|---|
| A | L12 exon 3, 79-bp two-cut deletion: detection limit and linearity | L12 cl6 (homozygous, 79 bp) | 0, 2, 5, 10, 25, 50, 75, 100 | 2 |
| B | Two near-identical deletions (79 vs 78 bp): allele identity | L12 cl6 + cl4 (78 bp) | 50/50; 25/25/50 WT; 45/5/50 WT | 2 |
| C | L07, single-base ssODN edit: intended-edit fraction | L07 H1.585 (homozygous edit, R2 1.000, no warnings) | 0, 2, 5, 10, 25, 50, 75, 100 | 2 |
| D | L08, two-variant ssODN edit | L08_5.744 (homozygous edit, R2 1.000) | 0, 5, 10, 25, 50, 100 | 2 |
| E | L09 exon 20, blocking-only conversion that other tools count as HDR | L09 Ex20 CI10 (98.4 % partial conversion, R2 0.999) | 0, 10, 25, 50, 100 | 2 |
| F1, F2 | Single-base indels: lowest detectable fraction | L03_5.1 (-1), L03_9.1 (+1 A) | 2, 5, 10, 25, 50 | 1 |
| WT | False-positive rate on parental DNA | parental line | 0 | 4 independent wells |

Source clones were chosen from the existing validation set: each was called homozygous by single-read analysis in the app with R2 of 0.97 or higher (one read direction only; a reverse read has not been taken). Confirm each source clone by out-out PCR (no larger second allele, no loss of heterozygosity) before it enters the series; a hidden second allele would move the truth.

## 3. Preparation
1. Amplify each locus from parental gDNA and from each source clone with the same primer pair (high-fidelity polymerase, 25 cycles; L12: L12_KO_Fw / L12_KO_Rev). Pool at least three reactions per template.
2. Purify on silica columns, elute in 10 mM Tris pH 8.0. Check one band on a gel or capillary system.
3. Quantify by fluorescence (dsDNA high-sensitivity dye), not absorbance, in triplicate. Convert to molarity: nM = ng/uL x 10^6 / (660 x amplicon bp). Dilute every template to 10 nM.
4. Prepare 10 % intermediate stocks (5 uL mutant + 45 uL WT) for each series that has levels below 10 %.
5. Build each mixture to 50 uL following `known-mixture-pipetting.csv` (no pipetting step below 2.5 uL). Two preparations per level are built from independent dilutions on different days.
6. Verify fractions independently where the amplicons differ in length (series A, B): capillary-electrophoresis peak areas, molarity-corrected. For series C, D, E, F amplicons are the same length: amplicon sequencing on a subset (levels 0, 5, 25, 50, 100, both preparations) gives the orthogonal check.

## 4. Sequencing
Submit 5 uL of each mixture (5 to 20 ng/uL) with the locus primer to the same facility, chemistry and instrument as routine clone QC (3730, BigDye v3). Forward reads and reverse reads of the same mixture go in adjacent wells of the same plate. Plate positions are randomised in `known-mixture-layout.csv`; each plate carries a water control in H12. Collect ab1 files with the well code in the file name.

## 5. Blinded analysis
The sealed key `known-mixture-key.csv` maps the code to the composition and is held by a person who does not run the analysis. The analyst gets the layout (code, direction, primer) and the parental control traces only. Analysis settings, tool versions and the scorecard thresholds are frozen and hashed (docs/scorecard.md) before the key is opened. All tools run on the same files, with the same guide/donor inputs.

## 6. What the series can and cannot show
- It measures accuracy, linearity, detection limit and false-positive rate on real traces with bench truth, in both read directions.
- Mixing after PCR removes PCR length bias by design; PCR bias is a separate effect (series A repeated with genomic DNA mixed before PCR is an optional extension) and is not part of the headline numbers.
- Sources are clone amplicons, so the allele spectrum is limited to the alleles listed above; it does not cover every indel.
- Five loci with two preparations of up to eight levels each do not support claims about loci not tested; generalisation relies on the in-silico set built from real controls of twelve loci.
