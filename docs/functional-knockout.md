# Functional knockout assessment

Assured QC reports two different numbers for a knockout, and they answer different questions.

- **Frameshift score** (also called the knockout score before 0.3.0): the share of alleles that carry an indel whose net length change is not a multiple of three, or that are 21 bp or longer. It looks at the trace only. It does not know where the edit lies in the gene, so it cannot tell a frameshift in the first exon from one in the last, and it counts a long in-frame deletion as a knockout whether or not it removes anything essential. It is a proxy.
- **Functional knockout score** (new, needs a gene record): the share of alleles predicted to lose function once the allele is placed in the coding sequence, translated, and judged on what happens to the transcript and the protein. It needs the gene record as a GenBank file.

Neither score is a measurement of protein. A knockout is confirmed by loss of protein (Western blot or targeted mass spectrometry) and, for knock-out designs that rely on transcript decay, by a fall in transcript downstream of the edit. The functional score is a sequence-based prediction of what that experiment will show, with the cases it cannot predict reported as uncertain instead of being scored as knockouts or as failures.

## Input

A GenBank record that contains the edited region and a coding sequence (CDS) feature with its exon structure, for example an NG_ (RefSeq genomic region) record from NCBI Nucleotide, or a gene-region export from Ensembl or SnapGene. The record is read in the browser or by the command line; nothing is fetched or uploaded.

- Command line: `--genbank gene.gb`, optionally `--transcript NM_000000` (otherwise the longest CDS is the primary isoform), `--regions "TIR domain:172-314,Death domain:67-122"` and `--epitope 200-250` (residue numbers of the protein in the record).
- Browser: choose the `.gb` file together with the other files (or in the options of the start page); regions and epitope are optional fields.
- API: `planDesign({ ..., genbank, codingOptions: { gene, transcript, regions, epitope } })`, then `analysePair` with `options.coding`.

- Design file: a design written by the Assured design app for a knockout or deletion carries the gene region and coding exons in its `coding` block (`docs/design-handoff.md`); no second upload is needed. A gene record given with the traces takes precedence.

A record that also carries a neighbouring gene is handled: only the coding sequences of the gene in which the reference window lies are modelled, and the others are named in a warning. A window that lies in no gene keeps every coding sequence.

The window that the traces cover (the design reference) is located in the record by k-mer seeds on either strand (identity of at least 0.9 required); if it is not found, or the record has no CDS, the functional score is omitted and the report says why, while the frameshift score is unchanged.

## What is computed for each allele

1. The allele sequence is compared to the reference window (common prefix and suffix) to find the edit, which is mapped to transcript coordinates (strand-aware, across exon junctions).
2. The edited transcript is rebuilt for each annotated CDS, translated, and the first stop codon is found.
3. Splice sites: if an edit removes or changes the canonical dinucleotides (GT at the 5' end, AG at the 3' end of an intron) at an exon boundary, two outcomes are evaluated, exon skipping and intron retention, and the allele is called only if both agree. A deletion that can be written in several places along a repeat is read where it leaves the canonical dinucleotide intact, because the allele sequence is the same and a GT or AG re-formed at the junction is a working site; such alleles are flagged as next to a splice site (the consensus around the site is changed, so splicing could still differ).
4. Classes and rules (all thresholds are configurable in `CODING_DEFAULTS`):

| Class | Meaning | Loss of function |
|---|---|---|
| `wt`, `no_coding_effect`, `synonymous` | the protein is unchanged | no |
| `lof_nmd` | premature stop predicted to trigger nonsense-mediated decay | yes |
| `lof_truncation` | premature stop that escapes NMD, truncation keeps less than 50% of the protein or lacks a region marked essential | yes |
| `truncation_uncertain` | stop that escapes NMD and keeps 50% or more of the protein | uncertain |
| `truncation_reinit_uncertain` | premature stop within 300 nt of the start with an in-frame ATG downstream that keeps at least 60% of the protein | uncertain |
| `lof_start_loss` / `start_loss_uncertain` | start codon destroyed; yes when no in-frame ATG downstream retains at least 60% | yes / uncertain |
| `lof_splice` / `splice_uncertain` | splice site destroyed; yes when both exon skipping and intron retention are loss of function | yes / uncertain |
| `inframe` | in-frame deletion, insertion or substitution | uncertain |
| `inframe_region_hit` | in-frame change that overlaps a region you marked essential | yes |
| `stop_loss`, `frameshift_no_stop` | extension past the normal stop | uncertain |
| `lof_deleted` | all coding exons are removed | yes |

5. NMD rules, from the exon-junction model (Nagy and Maquat 1998) as refined in human data (Lindeboom et al. 2016, 2019): a premature stop is predicted to trigger NMD unless it lies in the last coding exon, within 55 nt of the last exon-exon junction, within the first 150 coding nucleotides (start-proximal), in a coding exon longer than 407 nt, or the gene has a single coding exon. These rules are probabilistic: about a third of premature stops predicted to trigger NMD are not degraded in practice, so every NMD-based call carries the instruction to check the transcript and the protein.

## Clone calls

Alleles below 5% of the explained signal are not part of the call. The call is one of:

- `predicted_knockout`: every allele is predicted loss of function (NMD or clean loss), and nothing in the call is a known weak point.
- `predicted_knockout_with_caveats`: every allele is predicted loss of function, but at least one rests on a truncation that escapes NMD, a predicted splicing outcome, a lost start codon, an edit next to a splice site, a premature stop within 20 nt of the 55-nt boundary of the NMD rule, or an isoform that lacks the edited sequence. Each caveat is listed.
- `not_knockout`: wild type is present (5% or more of the signal) or an allele has no effect on the protein.
- `uncertain`: no allele is wild type, not all alleles are predicted loss of function (in-frame, small C-terminal truncations, ambiguous splicing).
- `unclear`: no allele above the floor.

The score is the share of the explained signal in each class (`functionalKoPct`, `functionalUncertainPct`), with the same profile-likelihood intervals as the other summaries.

In the decision layer a knockout is accepted only for `predicted_knockout`; every other functional call (`predicted_knockout_with_caveats`, `uncertain`, `not_knockout`, `unclear`) is held for a decision by hand. A clone is rejected only by the earlier wild-type rule (wild type is the major allele, 50% or more, or the genotype is homozygous wild type), not by the functional call. Every call lists what to run to confirm it: Western blot or mass spectrometry with an antibody or peptide outside the edited region, RT-qPCR downstream of the stop for NMD, N- and C-terminal antibodies for truncations, RT-PCR across the flanking exons for splice-site edits, and a functional assay for in-frame changes.

## Known limits

- This is a prediction. Residual protein after a frameshift knockout is common: about one third of verified knockout lines made by frameshifts still expressed protein from translation reinitiation or exon skipping in one proteomic survey (Smits et al. 2019), and aberrant mRNA or protein products arose in about half of a panel of edited lines (Tuladhar et al. 2019). The tool flags the mechanisms it can see in sequence (reinitiation sites, exon skipping, truncation), but cannot see exonic splicing enhancer disruption, internal ribosome entry or genetic compensation.
- The NMD rules are probabilistic (see above).
- In-frame alleles are always uncertain unless you mark essential regions; the tool does not predict whether a protein with a deleted stretch works.
- Cryptic splice sites are not predicted; only the canonical dinucleotides are checked. A deletion that spans an exon boundary is evaluated as skipping versus retention.
- Only the sequence in the record is modelled: a large deletion on the other allele that removes a primer site, or an off-target event, appears as apparent homozygosity and cannot be seen here.
- The thresholds (truncation 0.5 and 0.9, reinitiation 300 nt and 60%) are defaults chosen from the literature, not fitted to data; no laboratory protein measurements have been used to calibrate this score.

## Predicted protein products (figure)

Before a Western blot the question is which antibody can tell "no protein" from "a shorter protein". `src/proteinFigure.js` draws the predicted products of the alleles on the protein, from the allele table above, and shades the stretches where an epitope would see every product or the full-length protein only.

- A premature stop gives the residues kept before the first changed residue (solid), then the frameshifted tail up to the stop (hatched; dashed outline where NMD is expected or borderline, so the product may not be made), with its length and a calculated mass (about 0.11 kDa per residue; real gels run differently).
- An in-frame change gives the protein minus the removed residues. A splice-site loss gives its two outcomes: exon skipping (the residues encoded before the skipped exon are kept) and intron retention. Alleles that give the same product are one row.
- Alleles without a drawable product (start codon lost, no stop codon, all coding exons deleted) are listed as text. Wild-type and synonymous alleles draw nothing.
- Marked regions (`--regions`) and coding exons are drawn above the products. An epitope given with `--epitope` is drawn as a dashed box with the number of products that keep it.
- Orange stretch: an epitope inside it sees every product shown (so "no band" means nothing is made). Purple stretch: an epitope inside it sees the full-length protein only (a shorter product stays invisible). Residues are numbered in the protein of the gene record; vendors usually give UniProt numbers, which can differ by an offset when the isoforms differ.

Where it appears: under the allele table of each sample in the browser (SVG and 300 dpi PNG buttons) and in the HTML report, as `protein_products.svg` for the whole run (command line) and as a download after a run in the browser. It is a prediction from sequence, not a measurement.

## Tests

`test/coding.test.js` (record parser, plus- and minus-strand records, window location, stop positions against an independent translation, in-frame and last-exon rules, splice sites, start loss, isoforms, clone call), `test/coding-integration.test.js` (frameshift score unchanged without a record, agreement with an independent translation on the NR2F2 fixtures, decision, report and CSV) and `test/cli-gene.test.js` (command line end to end) and `test/protein-figure.test.js` (the figure of predicted products).

## References

- Nagy E, Maquat LE. A rule for termination-codon position within intron-containing genes: when nonsense affects RNA abundance. Trends Biochem Sci 1998;23:198-199.
- Lindeboom RGH, Supek F, Lehner B. The rules and impact of nonsense-mediated mRNA decay in human cancers. Nat Genet 2016;48:1112-1118.
- Lindeboom RGH, Vermeulen M, Lehner B, Supek F. The impact of nonsense-mediated mRNA decay on genetic disease, gene editing and cancer immunotherapy. Nat Genet 2019;51:1645-1651.
- Smits AH, Ziebell F, Joberty G, et al. Biological plasticity rescues target activity in CRISPR knock outs. Nat Methods 2019;16:1087-1093.
- Tuladhar R, Yeu Y, Piazza JT, et al. CRISPR-Cas9-based mutagenesis frequently provokes on-target mRNA misregulation. Nat Commun 2019;10:4056.
