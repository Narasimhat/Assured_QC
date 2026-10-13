# Which clones to carry forward

The decision (accept, hold, re-sequence, reject) says whether a result is reliable and what the genotype is. It cannot say which of several held clones to test first. `src/recommend.js` ranks the clones of a run by
how strongly the sequence predicts that no functional protein is left (the **knockout support**), and names a short panel to take to the bench. Nothing here changes a decision; a protein result does, through the
rule below.

## Knockout support

With a gene record, each allele of the functional assessment (`docs/functional-knockout.md`) gets a risk that it still makes a functional protein. The support of a clone is
100 x (1 - the share-weighted mean risk).

| Allele class | Risk |
| --- | --- |
| all coding exons deleted | 0.05 |
| premature stop, NMD predicted | 0.15 |
| the same, with a borderline NMD distance or next to a splice site | 0.30 |
| start codon lost | 0.30 |
| truncated protein that lacks the marked domain, NMD not expected | 0.40 |
| splice-site loss | the weaker of its two outcomes (exon skipping, intron retention) + 0.10 for a cryptic site |
| undecided prediction (in frame, no stop, reinitiation) | 0.70 |
| synonymous, no coding effect, short in-frame change | 0.95 |
| wild type | 1.00 |

The weights are ordinal. They put a clean NMD-predicted loss above a truncation that escapes decay above a call that rests on splicing; they are not measured probabilities, and a few points of difference
between two clones is not evidence. The best score a sequence alone can give is 85; a clone reaches 100 only when a protein result says the protein is absent. The tier follows the weakest allele: firm (risk up to
0.20), likely (0.35), possible (0.55), unsupported. Without a gene record the support is the frameshift proxy times the non-wild-type share, and the output says so. For SNP and tag edits it is the lower bound of
the intended-edit interval. Replace the weights with `options.risk` once `assured-qc calibrate` shows how often each class was borne out (it reports the absent-protein rate by allele class).

## Order and panel

Clones that were rejected or need re-sequencing are not ranked. The others are sorted by support, then by fewer distinct alleles (a homozygous clone is easier to read than a compound heterozygote), then by absence of
a trace-quality warning, then by confidence, then by name. Clones with the same genotype (the same alleles) form a group; the best of each group is a candidate, and the first three candidates (`--panel N`) are
**carry forward**. The rest are **reserve**; the other members of a group are **backup** (same genotype, so the same predicted outcome).

For a wanted edit (SNP, small tag, knock-in) the clones of one genotype are replicates, not backups: the clones carrying the edit are ranked one by one by the lower bound of the intended-edit interval, and the first `--panel N` of them are carry forward, so a plate with six clean homozygous clones can name all six. A clone whose sequence does not support the wanted edit (support tier unsupported, for example no intended edit, or a wild-type allele in a knockout) is **not recommended**; it is left out of the ranks and the panel. The reasons list the allele-level mechanism, and flags name what weakens
a ranking: a low-quality trace, a single read, low confidence in the edited share.

```
node cli/assured-qc.mjs run reads/ --design design.html --genbank gene.gb --regions "..." --out out/      # writes out/recommendation.csv, a section in plate.html, columns in results.csv
node cli/assured-qc.mjs recommend out/results.json [--observations protein.csv] [--panel 3] [--out recommendation.csv]
```

The browser shows the same ranking: a Priority column, the support and the reasons for the selected clone, and a download of the recommendation table.

## Run summary graph

`src/summaryFigure.js` draws a run on one graph, one row per sample in the order to act on them (ranked clones first): the decision (letter and colour, as in the plate view), what the clone is made of (a stacked bar of the alleles with their length changes and shares of the signal, and the signal no allele explains), the edited share (or the intended edit) with its 95% interval, the knockout support, and the plan. A `!` marks low trace quality or low confidence; hover text carries the full name and the reasons. It is shown above the table in the browser (SVG and 300 dpi PNG downloads), at the top of the HTML report, and written as `summary.svg` by `run`. Rows shrink for runs of more than 30 and 80 samples.

## Protein results

`--observations protein.csv` (the file that `calibrate` reads: `sample`, `protein` = absent / reduced / present / truncated, optional `rna`, `notes`) is applied to knockout and deletion samples only.

| Protein | Decision |
| --- | --- |
| absent | accept, if the sequence shows a clean knockout (edited at least 90%, wild type at most 10%); otherwise the decision stays and the reason says so |
| present | reject: the protein is still made, whatever the sequence predicts |
| truncated | hold: a shorter product is made; test whether it is functional |
| reduced | hold |

A clone whose protein is truncated or reduced is ranked below every unmeasured clone and marked **caution** (a protein is still made), and it is not counted toward the panel. The earlier decision is kept in `decision.before`, and the result is written into the reasons. Use an antibody whose epitope lies outside the edited region. An antibody that binds in the lost part reads `absent` for
a truncated product and proves nothing; two antibodies to different regions (one on each side of the edit) are better. The tool cannot check the epitope of the antibody you used.

## Limits

The ranking is a prediction from sequence, not a measurement. The weights have not been tested against protein data, and with few clones nothing can be tested. A clone with one read has an allele resolution that
rests on that read.


## KO score mode (no gene record)

Without a gene record the number in the support column is the KO score, as ICE and DECODR report it: the share of alleles with a frameshift or an indel of 21 bp or more, taken from the share of the signal that is not wild type. The summary graph, the plate view, the clone selection figure and the recommendation file name it "KO score" (with a gene record they say "Knockout support"). It says the change is likely to disrupt the gene, not that the protein is gone; a Western blot shows that.
Clones are ranked by it; ties (every clone with a large deletion between the cuts scores about 100) are broken by homozygous before compound heterozygous, then by a trace-quality warning, then by the confidence of the fit, then by the lower bound of the KO-score interval, then by name. Up to two clones with the same genotype count towards the panel, because independent clones with one deletion are a real replicate for a Western blot and a clean homozygous clone should come before a compound heterozygote with a noisy trace; the others are backups. With a gene record (functional prediction) one clone per genotype is taken, as before.
The default panel is three clones for a knockout and six for a wanted edit; `--panel N` sets it.


## HDR score mode (SNP, small tag)

For a wanted edit the number is the HDR score, which is what ICE reports as its HDR or KI score: the share of alleles that carry the intended edit (`summary.intendedEditPct`), with its 95% interval. Blocking changes without the edit are counted separately ("blocking-only") and are not part of it, which is where the app and ICE differ most: ICE counts any allele that carries donor bases as HDR.
The summary graph, the plate view, the clone selection figure, the report tables and the browser name it "HDR score". Clones are ranked by the lower bound of its interval, so a clone with a wide interval does not outrank a clean one on a point estimate.
A design with two ssODNs (one per guide, each with its own silent change) is read with both donors together (`docs/clone-selection.md`), so a clone made with either one scores as the intended edit. ICE takes one donor per run: a clone made with the other donor scores about 50 there, and a clone that has one donor's silent change on each allele (both alleles carry the edit) scores 70 to 80 with either, because each donor explains only part of the signal.
