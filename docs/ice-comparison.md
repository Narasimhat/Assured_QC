# Compare a run with ICE

After an analysis, expand **Compare with ICE**, then **Download ICE batch folder (ZIP)**. The package is created locally from the analysed trace pairs and design. It does not submit data to EditCo.

Extract the package. Each batch folder contains `ab1_files.zip` and `ice_definitions.xlsx`. In [ICE Batch Upload](https://ice.editco.bio/#/), select these two files. Keep the downloaded ICE results beside the batch. `sample_map.csv` and `sample_map.json` map each label to its original sample, read, donor model and Assured QC result; the JSON also includes input hashes and exclusions.

Two ssODNs produce two rows per read, retaining all analysed guides in both rows. They are alternative donor fits of the same sample, not a joint two-donor fit or independent replicates. Never sum or average their HDR scores. Identical donor sequences are exported once. Forward/reverse reads are separate rows with their matched controls (or the shared control used in Assured QC).

The five workbook headers and comma-separated multi-guide format match the [official ICE example package](https://ice.editco.bio/misc/ice_example_files.zip), checked 13 October 2026. Its batch template has one donor column and no nuclease column, so this exporter limits automatic preparation to explicit SpCas9 designs with 1–3 guides of 17–23 nt. Other nucleases need manual setup on ICE's sample-by-sample page. No guide sequence is invented for an inferred cut.

Donors include both arms, are limited to 16–300 DNA bases and need at least 15 aligned bases on each arm. ICE must validate the actual homology alignment. Unsupported or incomplete designs are listed as exclusions. Failed Assured QC results remain exportable for comparison with their original QC status in the map.

The default is 96 rows per batch; 20 is available for smaller uploads. All batches stay within 700 rows/files and a conservative 200 MB uncompressed AB1 limit. The outer ZIP is a folder download, not the AB1 archive ICE expects.

ICE donor HDR can differ from the intended amino-acid-change fraction, particularly where partial donor conversion carries the desired substitution. Compare fit quality, per-base signals and outcome definitions, not just headline percentages. Importing ICE result files for an automatic side-by-side comparison is a possible future extension.
