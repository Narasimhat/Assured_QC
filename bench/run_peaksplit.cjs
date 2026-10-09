
const fs = require("fs"), path = require("path");
const core = require(path.join(process.argv[2], "src", "core.js"));
const jobs = JSON.parse(fs.readFileSync(process.argv[3], "utf8")); const out = [];
const load = (p) => { const b = fs.readFileSync(p); return core.parseAB1(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength), path.basename(p)); };
const t0 = Date.now();
jobs.forEach((j, i) => {
  const t1 = Date.now(); let rec = { index: i };
  try {
    const control = load(j.control), edited = load(j.edited);
    // one guide per sample: the guide whose cut lies first in the control read
    let best = null;
    for (const g of j.guides) { try { const f = core.findGuide(control.seq, g, "SpCas9"); if (!best || f.cut < best.cut) best = { guide: g, cut: f.cut }; } catch (e) {} }
    if (!best) throw new Error("no guide found in control");
    const opts = j.kind === "KO" ? { maxDel: 100 } : {};
    const r = core.analyzeSample({ control, edited, guide: best.guide, nuclease: "SpCas9", donor: j.donor || "", opts });
    rec = { ...rec, ok: true, guideUsed: best.guide, wt: r.wtPct, indel: r.indelPct, ko: r.koScore, ki: r.kiScore, partialHDR: r.partialHDR, r2: r.r2, warnings: r.warnings,
      top: r.contribs.slice(0, 5).map((c) => [c.label, c.kind, Math.round(c.weight * 1000) / 1000]), donorSites: r.donorSites ? r.donorSites.map((d) => ({ pos: d.pos, edited: d.edited })) : null };
  } catch (e) { rec = { ...rec, ok: false, error: String(e.message || e).slice(0, 200) }; }
  rec.seconds = (Date.now() - t1) / 1000; out.push(rec);
  if ((i + 1) % 20 === 0) console.error(`${i + 1}/${jobs.length} ${Math.round((Date.now() - t0) / 1000)}s`);
});
fs.writeFileSync(process.argv[4], JSON.stringify(out));
