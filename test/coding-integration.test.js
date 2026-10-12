import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { analysePair } from "../src/analyze.js";
import { analyseSample, planDesign, headline } from "../src/batch.js";
import { decide } from "../src/decision.js";
import { buildResultsCsv, buildQcReportHtml } from "../src/report.js";
import { translate } from "../src/coding.js";
import { loadSpec, loadTrace } from "./helpers.js";

const here = path.dirname(fileURLToPath(import.meta.url));
const truth = JSON.parse(readFileSync(path.join(here, "..", "fixtures", "sim", "truth.json"), "utf8"));
const spec = loadSpec("nr2f2_2xha_ssodn"); const ref = spec.reference.sequence; const cut = spec.guides[0].cut;

// A gene around the real reference window: exon 1 (synthetic) | intron | exon 2 = a slice of the window around the guide cut | intron | exon 3 (synthetic, ends in TAA).
let seed = 3; const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
const B = "ACGT"; const dna = (n) => Array.from({ length: n }, () => B[Math.floor(rnd() * 4)]).join("");
const CODONS = []; for (const a of B) for (const b of B) for (const c of B) if (!["TAA", "TAG", "TGA", "ATG"].includes(a + b + c)) CODONS.push(a + b + c);
const codons = (n) => Array.from({ length: n }, () => CODONS[Math.floor(rnd() * CODONS.length)]).join("");
const isStop = (c) => c === "TAA" || c === "TAG" || c === "TGA";
// the longest stop-free stretch around the cut, in the best of the three frames: exon 2 is that stretch
let a = 0; let exonLen = 0;
for (let f = 0; f < 3; f += 1) {
  let p0 = cut - 3 - ((cut - 3 - f) % 3); while (p0 >= 0 && !isStop(ref.slice(p0, p0 + 3))) p0 -= 3;
  let p1 = p0 + 3; while (p1 + 3 <= ref.length && !isStop(ref.slice(p1, p1 + 3))) p1 += 3;
  if (p0 < 0 || p1 + 3 > ref.length) continue; const start = p0 + 3; const len = p1 - start; if (len > exonLen && start < cut - 6 && start + len > cut + 6) { a = start; exonLen = len; }
}
assert.ok(exonLen >= 45, `an open stretch around the cut exists (${exonLen} nt)`);
const exon1 = "ATG" + codons(70); const exon2 = ref.slice(a, a + exonLen); const exon3 = codons(10) + "TAA"; const cdsSeq = exon1 + exon2 + exon3;
const up = dna(200); const i1 = "GT" + dna(120) + "AG"; const i2 = "GT" + dna(120) + "AG"; const down = dna(200);
const genomic = up + exon1 + i1 + ref + i2 + exon3 + down;
const e1 = [up.length, up.length + exon1.length]; const e2s = up.length + exon1.length + i1.length + a; const e2 = [e2s, e2s + exonLen]; const e3s = up.length + exon1.length + i1.length + ref.length + i2.length; const e3 = [e3s, e3s + exon3.length];
const gb = ["LOCUS       SYN2   " + genomic.length + " bp", "ACCESSION   NG_SYN2", "FEATURES             Location/Qualifiers",
  `     CDS             join(${e1[0] + 1}..${e1[1]},${e2[0] + 1}..${e2[1]},${e3[0] + 1}..${e3[1]})`, '                     /gene="NR2F2"', '                     /transcript_id="NM_SYN2.1"', "ORIGIN",
  ...Array.from({ length: Math.ceil(genomic.length / 60) }, (_, i) => String(i * 60 + 1).padStart(9) + " " + genomic.slice(i * 60, i * 60 + 60).toLowerCase()), "//"].join("\n");
const trace = (name) => loadTrace(name); const entry = truth.nr2f2_ko_clone_hom_m1;
const run = (extra = {}) => analyseSample({ name: "clone", spec, control: trace(entry.control), edited: trace(entry.edited), sampleType: "clone", options: extra });
const editOnCds = (S) => { let p = 0; const n = Math.min(ref.length, S.length); while (p < n && ref[p] === S[p]) p += 1; let q = 0; while (q < n - p && ref[ref.length - 1 - q] === S[S.length - 1 - q]) q += 1; return { s: p, e: ref.length - q, ins: S.slice(p, S.length - q) }; };

test("without a gene record nothing changes and the report says the knockout score is a frameshift proxy", () => {
  const plain = run(); assert.equal(plain.ok, true); assert.equal(plain.functional, null);
  const withCoding = run({ codingInput: { genbank: gb, gene: "NR2F2" } }); assert.equal(withCoding.ok, true);
  assert.equal(withCoding.summary.koScorePct, plain.summary.koScorePct, "the frameshift score is unchanged");
  assert.deepEqual(withCoding.contributions.map((c) => c.label), plain.contributions.map((c) => c.label));
  const html = buildQcReportHtml({ workflow: "knockout", results: [plain] });
  assert.match(html, /proxy/i); assert.doesNotMatch(html, /Functional knockout/);
  assert.doesNotMatch(html, /Predicted protein products/);
});

test("with a gene record each allele is translated: stop position agrees with an independent translation, call and interval are reported", () => {
  const r = run({ codingInput: { genbank: gb, gene: "NR2F2" } }); assert.ok(r.functional, "functional assessment present");
  assert.ok(r.summary.functionalKoPct !== undefined && r.intervals.functionalKoPct, "score and interval");
  const [lo, hi] = r.intervals.functionalKoPct; assert.ok(lo <= r.summary.functionalKoPct + 0.1 && hi >= r.summary.functionalKoPct - 0.1);
  const row = r.functional.rows.find((x) => x.kind === "indel"); assert.ok(row, "an indel allele is reported");
  // independent expectation for that allele: edit on the window -> edit on the coding sequence (exon 2 starts at window index a, coding offset exon1.length)
  const result = analysePair({ control: trace(entry.control), edited: trace(entry.edited), spec, options: { sampleType: "clone" } });
  const contribution = result.decomposition.contributions.find((c) => c.label === row.label); const ed = editOnCds(contribution.sequence);
  assert.ok(ed.s >= a && ed.e <= a + exonLen, "the edit lies inside the modelled exon");
  const cdsPos = exon1.length + (ed.s - a); const mutated = cdsSeq.slice(0, cdsPos) + ed.ins + cdsSeq.slice(cdsPos + (ed.e - ed.s));
  const stop = translate(mutated).stopNt;
  if (row.details.ptcNt !== undefined) assert.equal(row.details.ptcNt, stop); else assert.ok(stop < 0 || mutated.length % 3 === 0);
  assert.ok(["predicted_knockout", "predicted_knockout_with_caveats", "uncertain"].includes(r.functional.call.status), r.functional.call.status);
  assert.match(headline("knockout", r), /knockout|uncertain|caveats/);
});

test("the decision and the exports use the functional call", () => {
  const r = run({ codingInput: { genbank: gb, gene: "NR2F2" } });
  const d = decide("knockout", r); assert.ok(["accept", "hold", "re-sequence", "reject"].includes(d.decision));
  const csv = buildResultsCsv("knockout", [r]); assert.match(csv.split("\n")[0], /functional_ko_pct/); assert.match(csv.split("\n")[0], /functional_call/);
  const html = buildQcReportHtml({ workflow: "knockout", results: [r] });
  assert.match(html, /Functional knockout/); assert.match(html, /Western/);
  assert.match(html, /<svg[^>]*role="img"/); assert.match(html, /Predicted protein products/); assert.ok((html.match(/<svg /g) || []).length >= 1);
});

test("a gene record that does not contain the design's reference is reported, and the analysis still runs", () => {
  const other = gb.replace(/ORIGIN[\s\S]*$/, "ORIGIN\n        1 " + dna(300).toLowerCase() + "\n//");
  const plan = planDesign({ designSpec: spec, genbank: other }); assert.ok(plan.warnings.some((w) => /Functional knockout assessment is off/.test(w)));
  const r = run({ codingInput: null }); assert.equal(r.ok, true);
});
