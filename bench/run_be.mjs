import { readFileSync, writeFileSync } from "node:fs";
import { runJobs } from "../src/runJobs.js";
const jobs = JSON.parse(readFileSync(process.argv[2], "utf8")); const out = [];
const file = (p, name) => { const b = readFileSync(p); return { name, buffer: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) }; };
jobs.forEach((j, i) => {
  const r = runJobs({ workflow: "base_edit", manual: { guides: j.guides.join(", "), gene: j.project }, baseEditor: { editor: j.editor, window: [4, 8], target: j.target }, sampleType: "pool", options: { baseEdit: process.env.QC_BE_OPTIONS ? JSON.parse(process.env.QC_BE_OPTIONS) : {} }, control: file(j.control, "c.ab1"), samples: [file(j.edited, "s.ab1")] });
  const s = r.ok ? r.results[0] : null;
  out.push(r.ok && s.ok ? { i, ok: true, edits: s.baseEdits.positions, summary: s.summary } : { i, ok: false, error: r.error || (s && s.error) });
});
writeFileSync(process.argv[3], JSON.stringify(out));
