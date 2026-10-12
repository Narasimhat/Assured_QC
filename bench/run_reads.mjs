// Runs each job three ways: forward read alone, reverse read alone, and both fitted together. Usage: node bench/run_reads.mjs jobs.json out.json [first last]
import { readFileSync, writeFileSync } from "node:fs";
import { parseAbif } from "../src/abif.js";
import { planDesign } from "../src/batch.js";
import { analyseReads } from "../src/analyze.js";

const [jobsPath, outPath, first = "0", last = "1000000"] = process.argv.slice(2);
const jobs = JSON.parse(readFileSync(jobsPath, "utf8")).slice(Number(first), Number(last));
const trace = (p) => { const b = readFileSync(p); return parseAbif(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength), p); };
const brief = (r) => (r.ok ? { ok: true, edited: r.summary.editedPct, intended: r.summary.intendedEditPct, ko: r.summary.koScorePct, interval: r.intervals?.editedPct, r2: r.decomposition.r2, warnings: r.warnings.length, reads: r.reads || null, agreement: r.agreement || null } : { ok: false, error: r.error });
const out = [];
jobs.forEach((job, index) => {
  const row = { index, label: job.label };
  try {
    const cF = trace(job.control), eF = trace(job.edited), cR = trace(job.control_rev), eR = trace(job.edited_rev);
    const plan = planDesign({ designSpec: null, control: cF, guides: job.guides.join(", "), donor: job.donor || "", gene: job.project });
    const options = { sampleType: job.sample_type || "clone" };
    const t0 = Date.now();
    row.fwd = brief(analyseReads({ reads: [{ name: "F", control: cF, edited: eF }], spec: plan.spec, options }));
    row.rev = brief(analyseReads({ reads: [{ name: "R", control: cR, edited: eR }], spec: plan.spec, options }));
    row.joint = brief(analyseReads({ reads: [{ name: "F", control: cF, edited: eF }, { name: "R", control: cR, edited: eR }], spec: plan.spec, options }));
    row.seconds = (Date.now() - t0) / 1000;
  } catch (error) { row.error = String(error.stack || error).slice(0, 300); }
  out.push(row);
});
writeFileSync(outPath, JSON.stringify(out));
console.log(`done ${out.length}`);
