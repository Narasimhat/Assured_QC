# Calibrating the functional knockout call

The functional knockout call (`docs/functional-knockout.md`) is a prediction from sequence. Its thresholds (NMD rules, truncation 50% and 90%, reinitiation 300 nt and 60%) are
literature defaults and have not been fitted to any protein data. This is how to measure how often a call is borne out.

1. Run the clones with a gene record (`--genbank`, or a design file with a coding block); keep `results.json` from the output folder.
2. For each clone with a Western blot or mass-spectrometry result, write one row (`docs/calibration-template.csv`): `sample` (the name used in the run), `protein`
   (`absent`, `reduced`, `present` or `truncated`), and optionally `rna` (`decayed` or `stable`, from RT-qPCR downstream of the stop) and `splicing` (free text from RT-PCR). Use an antibody
   whose epitope lies outside the edited region, with the parental line as positive control.
3. `node cli/assured-qc.mjs calibrate results.json observations.csv --out calibration.md`.

The report gives, for each call, how many clones had no protein, with Wilson 95% intervals; the sensitivity and specificity of the knockout call; and the same by the most abundant allele class
(NMD, truncation, splice site, in-frame), which shows which rule fails. Clones with a measured result that are not in the run, or protein words the tool does not recognise, are listed, not guessed.

The report does not change any threshold. With fewer than about 20 clones the intervals are wide and the report says so; a threshold should only be moved when a class has enough clones to show a rule is wrong,
and the change should then be tested on clones that were not used to choose it.

The same table (`sample`, `protein`, optional `rna`) can be given to `run` or `recommend` with `--observations`: a clone whose protein is absent is then accepted, one whose protein is present is rejected (`docs/recommendation.md`).
