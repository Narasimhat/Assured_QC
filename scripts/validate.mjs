// Compares this engine, the known allele mixtures and (where available) Synthego ICE on the synthetic traces in fixtures/sim.
// Run: node scripts/validate.mjs [out.csv]. Prints a table and the mean absolute error against the known mixture.
import { readFileSync, writeFileSync } from "node:fs";
import path from "node:path";
import { parseAbif } from "../src/abif.js";
import { analysePair } from "../src/analyze.js";
import { loadSpec, loadTruth, simDir } from "../test/helpers.js";

const truth = loadTruth();
const oracle = JSON.parse(readFileSync(path.join(simDir, "oracle_ice.json"), "utf8")).results;
const cache = new Map();
const trace = (name) => { if (!cache.has(name)) cache.set(name, parseAbif(readFileSync(path.join(simDir, name)), name)); return cache.get(name); };
const specFor = (name) => loadSpec(name.startsWith("nr2f2") ? "nr2f2_2xha_ssodn" : "apoe_r176c_snp");
const sizeOf = (label, entry) => { if (label === "del_between_cuts") return -(entry.cuts.g2 - entry.cuts.g1); const m = /(del|ins)([+-]?)(\d+)/.exec(label); return m ? (m[1] === "del" ? -1 : 1) * Number(m[3]) : null; };
const isKo = (size) => size !== null && (size % 3 !== 0 || Math.abs(size) >= 21);
const pct = (value) => Math.round(value * 10) / 10;

const rows = [];
for (const [name, entry] of Object.entries(truth)) {
  const result = analysePair({ control: trace(entry.control), edited: trace(entry.edited), spec: specFor(name), options: { sampleType: /clone/.test(entry.kind) || /het|hom|compound|null/.test(name) ? "clone" : "pool" } });
  if (!result.ok) { rows.push({ sample: name, error: result.error }); continue; }
  const sum = (test) => 100 * entry.alleles.filter(test).reduce((total, a) => total + a.fraction, 0);
  const tWt = sum((a) => a.label === "wt"), tIntended = sum((a) => a.label === "hdr"), tBlock = sum((a) => a.label === "silent_only");
  const tKo = sum((a) => isKo(sizeOf(a.label, entry)));
  const ice = oracle[name] && oracle[name].summary;
  const s = result.summary;
  rows.push({
    sample: name, truth_wt: pct(tWt), own_wt: s.wtPct, truth_edited: pct(100 - tWt), own_edited: s.editedPct, ice_edited: ice ? ice.ice : "",
    truth_intended: pct(tIntended), own_intended: s.intendedEditPct, ice_hdr: ice ? pct(ice.hdr_pct) : "", truth_blocking_only: pct(tBlock), own_blocking_only: s.partialConversionPct,
    truth_ko: pct(tKo), own_ko: s.koScorePct, ice_ko: ice ? ice.ko_score : "", own_r2: result.decomposition.r2 === null ? "" : Math.round(result.decomposition.r2 * 1000) / 1000,
  });
}
const cols = Object.keys(rows.find((r) => !r.error));
const csv = [cols.join(","), ...rows.filter((r) => !r.error).map((r) => cols.map((c) => r[c]).join(","))].join("\n") + "\n";
if (process.argv[2]) writeFileSync(process.argv[2], csv);
console.log(csv);
const mae = (own, ref, filter = () => true) => { const used = rows.filter((r) => !r.error && filter(r) && r[ref] !== "" && r[own] !== ""); return { n: used.length, mae: pct(used.reduce((t, r) => t + Math.abs(r[own] - r[ref]), 0) / Math.max(1, used.length)), max: pct(Math.max(0, ...used.map((r) => Math.abs(r[own] - r[ref])))) }; };
const bothHaveIce = (r) => r.ice_edited !== "";
console.log(JSON.stringify({
  edited_own_vs_truth: mae("own_edited", "truth_edited"), edited_own_vs_truth_on_ice_samples: mae("own_edited", "truth_edited", bothHaveIce), edited_ice_vs_truth_on_ice_samples: mae("ice_edited", "truth_edited", bothHaveIce),
  intended_own_vs_truth: mae("own_intended", "truth_intended"), intended_own_vs_truth_on_ice_samples: mae("own_intended", "truth_intended", bothHaveIce), hdr_ice_vs_truth: mae("ice_hdr", "truth_intended", bothHaveIce),
  ko_own_vs_truth: mae("own_ko", "truth_ko"), ko_own_vs_truth_on_ice_samples: mae("own_ko", "truth_ko", bothHaveIce), ko_ice_vs_truth: mae("ice_ko", "truth_ko", bothHaveIce),
  errors: rows.filter((r) => r.error),
}, null, 1));
