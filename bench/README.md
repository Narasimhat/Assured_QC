# Benchmark harness

One job file, one command, every comparator that can be scripted, one result shape.

    python bench/run_all.py jobs.json out/ --tools assured_qc,ice,peaksplit,traceedit
    python bench/score.py jobs.json out/            # scorecard (added with the truth data)

Paths to the comparators are set in `bench/config.json` (copy `bench/config.example.json`; the real file is git-ignored) (node, ICE checkout + its python, peaksplit checkout, TraceEdit checkout).

## Job format
`{project, label, control, edited, guides: [spacer, ...], donor: "" | ssODN, sample_type: "clone" | "pool", kind: "KO" | "SNP", truth: {...}}`
`control` and `edited` are .ab1 paths. `truth` is optional and used only by the scorer.

## Result format (per job, per tool)
`{ok, error, edited_pct, intended_pct, ko_pct, alleles: [{label, size, pct}], r2, seconds}`

## Comparator status on this host (Windows, no Linux compute attached)
| Tool | Scripted here | Notes |
|---|---|---|
| Assured QC | yes | node, 0.4 s/pair |
| ICE 1.2.0 | yes | local checkout, 22 s/pair on this host (parallel slices); non-commercial licence, used as a comparator only |
| peaksplit | yes | local checkout, one guide per run |
| Assured_TraceEdit | yes | local checkout; refuses pairs failing its early-read rule (5/29 in the smoke test) |
| Tracy | no | the bioconda package has no Windows build; needs a Linux compute target |
| DECODR, TIDE/TIDER, SeqScreener, EditR | no | web-only interfaces with no documented batch API; run by hand on a subset if a figure is needed |

## Public benchmark data
No public set of raw .ab1 traces with known allele fractions was found in the literature searched (Aoki et al. 2024 used artificial templates built in the lab and report values in figures and tables; no raw-trace deposit located). Truth therefore comes from (1) the in-silico generator in `bench/simulate.js` built from real control peaks, (2) the known-mixture protocol in `docs/known-mixture-protocol.md` (wet-lab, not yet run), and (3) an NGS subset of real clones, if the user can supply one.

## Scoring and secondary tables

    python bench/score.py jobs.json out/ prefix --split test      # frozen scorer (hash in FROZEN.json)
    python bench/extras.py jobs.json extras_in/ prefix --split test # intervals, confidence tiers, decisions (needs aq_out.json and dec_out.json from bench/run_decisions.mjs)

`bench/ENGINE_FREEZE.json` records the engine file hashes used for the held-out run in `docs/results.md`.
