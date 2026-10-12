import { readFileSync, writeFileSync } from "node:fs";
import { runJobs } from "../src/runJobs.js";
import { decide } from "../src/decision.js";
const [jobsPath, outPath] = process.argv.slice(2); const jobs = JSON.parse(readFileSync(jobsPath, "utf8")); const out = [];
const file = (p, name) => { const b = readFileSync(p); return { name, buffer: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) }; };
jobs.forEach((job, index) => {
  const workflow = job.donor ? "snp" : (job.guides.length > 1 ? "deletion" : "knockout");
  const r = runJobs({ workflow, manual: { guides: job.guides.join(", "), donor: job.donor || "", gene: job.project }, sampleType: job.sample_type || "clone", control: file(job.control, "c.ab1"), samples: [file(job.edited, "s.ab1")] });
  const s = r.ok ? r.results[0] : null;
  if (!r.ok || !s.ok) { out.push({ index, ok: false, error: r.error || s.error, decision: decide(workflow, { ok: false, error: r.error || (s && s.error) }) }); return; }
  out.push({ index, ok: true, decision: decide(workflow, s), category: s.genotype?.category || null, tier: s.confidence?.tier });
});
writeFileSync(outPath, JSON.stringify(out));
