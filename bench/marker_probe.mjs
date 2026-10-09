
import { readFileSync, writeFileSync } from "node:fs";
import { parseAbif } from "../src/abif.js";
import { planDesign } from "../src/batch.js";
import { analysePair } from "../src/analyze.js";
import { markerGain, correctedShare } from "../src/markerHeights.js";
const [jobsPath, outPath, only] = process.argv.slice(2); const jobs = JSON.parse(readFileSync(jobsPath, "utf8")); const out = [];
const load = (p) => { const b = readFileSync(p); return parseAbif(b, p.split(/[\\/]/).pop()); };
jobs.forEach((job, index) => {
  if (!job.donor || (job.sample_type || "clone") !== "clone") return;
  if (only && !new RegExp(only).test(job.project)) return;
  try {
    const control = load(job.control); const edited = load(job.edited);
    const plan = planDesign({ control, guides: job.guides.join(", "), donor: job.donor, gene: job.project, editedForCut: [edited] }); if (plan.error) return;
    const r = analysePair({ control, edited, spec: plan.spec, options: { sampleType: "clone" } }); if (!r.ok) return;
    const rows = r.markerReadout.filter((m) => m.covered).map((m) => { const g = markerGain(r.prepared.control, m.readIndex, m.marker.ref, m.marker.alt); return { pos: m.marker.pos, change: m.marker.change, raw: m.altNet, g, corrected: correctedShare(m.altNet, g) }; });
    out.push({ index, project: job.project, label: job.label, markers: rows, intended: r.summary.intendedEditPct, edited: r.summary.editedPct, partial: r.summary.partialConversionPct });
  } catch (e) { /* skip */ }
});
writeFileSync(outPath, JSON.stringify(out)); console.log(out.length);
