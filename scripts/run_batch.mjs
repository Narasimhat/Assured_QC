// Runs one batch exactly as the app does (same engine, same report builder) from the command line.
// Usage: node scripts/run_batch.mjs config.json
// config: { title, workflow, design, guides: "spacer1, spacer2", donor: "", gene, sampleType: "clone"|"pool", control: path, samples: [{name, path}], out: "folder" }
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import path from "node:path";
import { runJobs } from "../src/runJobs.js";
import { buildQcReportHtml, buildResultsCsv } from "../src/report.js";

const cfg = JSON.parse(readFileSync(process.argv[2], "utf8"));
const file = (p, name) => { const bytes = readFileSync(p); return { name, buffer: bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) }; };
const out = runJobs({
  workflow: cfg.workflow, manual: { guides: cfg.guides, donor: cfg.donor || "", gene: cfg.gene || "" }, sampleType: cfg.sampleType || "clone",
  control: file(cfg.control, path.basename(cfg.control)), samples: cfg.samples.map((s) => file(s.path, s.name)),
});
mkdirSync(cfg.out, { recursive: true });
if (!out.ok) { console.error(out.error); process.exit(1); }
writeFileSync(path.join(cfg.out, "results_full.json"), JSON.stringify(out));
writeFileSync(path.join(cfg.out, "report.html"), buildQcReportHtml({ title: cfg.title, workflow: cfg.workflow, design: cfg.design, results: out.results, notes: out.notes, generatedAt: new Date() }));
writeFileSync(path.join(cfg.out, "results.csv"), buildResultsCsv(cfg.workflow, out.results));
console.log(JSON.stringify(out.results.map((r) => ({ name: r.name, ok: r.ok, error: r.error, orientation: r.orientation, summary: r.summary, genotype: r.genotype?.summary, warnings: r.warnings?.length }))));
