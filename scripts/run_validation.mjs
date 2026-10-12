// Runs the QC engine on a list of control/edited pairs described in a JSON file, for comparison with another tool (ICE).
// Usage: node scripts/run_validation.mjs jobs.json results.json
// jobs.json: [{ project, label, control, edited, guides: [spacer, ...], donor: "" | "ssODN sequence", sample_type: "clone" | "pool" }]
// (control and edited are file paths). Guides and donor are entered the way a user would type them in the app.
import { readFileSync, writeFileSync } from "node:fs";
import { runJobs } from "../src/runJobs.js";

const [jobsPath, outPath] = process.argv.slice(2);
const jobs = JSON.parse(readFileSync(jobsPath, "utf8"));
const file = (path, name) => { const bytes = readFileSync(path); return { name, buffer: bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) }; };
const results = [];
for (const [index, job] of jobs.entries()) {
  const started = Date.now();
  const workflow = job.donor ? "snp" : (job.guides.length > 1 ? "deletion" : "knockout");
  let record;
  try {
    const out = runJobs({
      workflow, manual: { guides: job.guides.join(", "), donor: job.donor || "", gene: job.project }, sampleType: job.sample_type || "clone",
      control: file(job.control, "control.ab1"), samples: [file(job.edited, "sample.ab1")], options: process.env.QC_OPTIONS ? JSON.parse(process.env.QC_OPTIONS) : {},
    });
    if (!out.ok) record = { ok: false, error: out.error };
    else {
      const r = out.results[0];
      if (!r.ok) record = { ok: false, error: r.error, mixedSignal: /mixed/.test(r.error || "") };
      else record = {
        ok: true, workflow, planNotes: out.notes, orientation: r.orientation, summary: r.summary, intervals: r.intervals, quality: r.quality,
        genotype: r.genotype ? { category: r.genotype.category, summary: r.genotype.summary, flags: r.genotype.flags, alleles: r.genotype.alleles } : null,
        contributions: r.contributions.slice(0, 8), warnings: r.warnings,
      };
    }
  } catch (error) { record = { ok: false, error: `exception: ${error.message}` }; }
  record.index = index; record.seconds = (Date.now() - started) / 1000;
  results.push(record);
  if ((index + 1) % 20 === 0) console.log(`${index + 1}/${jobs.length}`);
}
writeFileSync(outPath, JSON.stringify(results));
console.log(`done: ${results.filter((r) => r.ok).length} analysed, ${results.filter((r) => !r.ok).length} not analysed`);
