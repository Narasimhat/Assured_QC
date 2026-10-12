import { readFileSync, writeFileSync } from "node:fs";
import { runJobs } from "../src/runJobs.js";
const [jobsPath, outPath] = process.argv.slice(2); const jobs = JSON.parse(readFileSync(jobsPath, "utf8")).filter((j) => j.donor);
const file = (p, name) => { const b = readFileSync(p); return { name, buffer: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) }; };
const rows = [];
for (const job of jobs) {
  try {
    const out = runJobs({ workflow: "snp", manual: { guides: job.guides.join(", "), donor: job.donor, gene: job.project }, sampleType: "clone", control: file(job.control, "control.ab1"), samples: [file(job.edited, "sample.ab1")], options: {} });
    const r = out.ok ? out.results[0] : null; if (!r || !r.ok) { rows.push({ set: job.project, label: job.label, ok: false }); continue; }
    rows.push({ set: job.project, label: job.label, ok: true, intendedPct: r.summary.intendedEditPct, partialPct: r.summary.partialConversionPct, wtPct: r.summary.wtPct, markers: (r.markers || []).map((m) => ({ pos: m.position, change: m.change, role: m.role, covered: m.covered, alt: m.altEdited, ctl: m.altControl, net: m.altNet, ref: m.refEdited })) });
  } catch (e) { rows.push({ set: job.project, label: job.label, ok: false, err: e.message }); }
}
writeFileSync(outPath, JSON.stringify(rows)); console.log(rows.filter((r) => r.ok).length, "of", rows.length);