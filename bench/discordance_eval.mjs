import { readFileSync, writeFileSync } from "node:fs";
import { runJobs } from "../src/runJobs.js";
const [jobsPath, outPath] = process.argv.slice(2); const jobs = JSON.parse(readFileSync(jobsPath, "utf8"));
const file = (p, name) => { const b = readFileSync(p); return { name, buffer: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) }; };
const rows = [];
for (const job of jobs) {
  const workflow = job.donor ? "snp" : (job.guides.length > 1 ? "deletion" : "knockout"); let rec = { set: job.project, label: job.label, workflow };
  try {
    const out = runJobs({ workflow, manual: { guides: job.guides.join(", "), donor: job.donor || "", gene: job.project }, sampleType: "clone", control: file(job.control, "control.ab1"), samples: [file(job.edited, "sample.ab1")], options: {} });
    const r = out.ok ? out.results[0] : null;
    if (!r || !r.ok) rec.ok = false; else rec = { ...rec, ok: true, editedPct: r.summary.editedPct, intendedPct: r.summary.intendedEditPct, koPct: r.summary.koScorePct, before: r.quality.discordanceBefore, nearCut: r.quality.discordanceNearCut ?? null, after: r.quality.discordanceAfter, r2: r.quality.r2 ?? null, upstreamRms: r.quality.upstreamRms, warned: r.warnings.some((w) => /already differ/.test(w)), nWarnings: r.warnings.length, conf: r.confidence && r.confidence.tier };
  } catch (e) { rec.ok = false; rec.err = e.message; }
  rows.push(rec);
}
writeFileSync(outPath, JSON.stringify(rows)); console.log(rows.filter((r) => r.ok).length, "of", rows.length);