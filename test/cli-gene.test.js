import test from "node:test";
import assert from "node:assert/strict";
import { mkdtempSync, copyFileSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { main as cli } from "../cli/assured-qc.mjs";
import { translate } from "../src/coding.js";

const here = path.dirname(fileURLToPath(import.meta.url)); const root = path.join(here, "..");
const spec = JSON.parse(readFileSync(path.join(root, "fixtures", "design", "nr2f2_2xha_ssodn.design.json"), "utf8")); const ref = spec.reference.sequence; const cut = spec.guides[0].cut;
let seed = 5; const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
const B = "ACGT"; const dna = (n) => Array.from({ length: n }, () => B[Math.floor(rnd() * 4)]).join("");
const CODONS = []; for (const a of B) for (const b of B) for (const c of B) if (!["TAA", "TAG", "TGA", "ATG"].includes(a + b + c)) CODONS.push(a + b + c);
const codons = (n) => Array.from({ length: n }, () => CODONS[Math.floor(rnd() * CODONS.length)]).join("");
const isStop = (c) => c === "TAA" || c === "TAG" || c === "TGA";
let a = 0; let len = 0;
for (let f = 0; f < 3; f += 1) { let p0 = cut - 3 - ((cut - 3 - f) % 3); while (p0 >= 0 && !isStop(ref.slice(p0, p0 + 3))) p0 -= 3; let p1 = p0 + 3; while (p1 + 3 <= ref.length && !isStop(ref.slice(p1, p1 + 3))) p1 += 3; if (p0 >= 0 && p1 - p0 - 3 > len && p0 + 3 < cut - 6 && p1 > cut + 6) { a = p0 + 3; len = p1 - p0 - 3; } }
const exon1 = "ATG" + codons(70);
const a2 = a + 3 * Math.floor(Math.max(0, cut - a - 90) / 3); const L2 = Math.floor((a + len - a2) / 3) * 3; const exon2 = ref.slice(a2, a2 + L2);
// a middle exon of random codons, drawn until a 1-nt deletion at the cut puts the stop more than 60 nt upstream of its junction (NMD expected)
const k = cut - a2; const e2del = exon2.slice(0, k - 1) + exon2.slice(k); let exon3 = ""; let tries = 0;
for (; tries < 200; tries += 1) { exon3 = codons(70); const stop = translate(exon1 + e2del + exon3).stopNt; if (stop > exon1.length + e2del.length && stop < exon1.length + e2del.length + exon3.length - 60) break; }
assert.ok(tries < 200, "found a middle exon");
const exon4 = codons(10) + "TAA"; const up = dna(200); const i1 = "GT" + dna(120) + "AG"; const i2 = "GT" + dna(120) + "AG"; const i3 = "GT" + dna(120) + "AG";
const pre = ref.slice(0, a2); const post = ref.slice(a2 + L2);
const genomic = up + exon1 + i1 + pre + exon2 + post + i2 + exon3 + i3 + exon4 + dna(200);
const e2s = up.length + exon1.length + i1.length + a2; const e3s = up.length + exon1.length + i1.length + ref.length + i2.length; const e4s = e3s + exon3.length + i3.length;
const gb = ["LOCUS       SYN3   " + genomic.length + " bp", "ACCESSION   NG_SYN3", "FEATURES             Location/Qualifiers", `     CDS             join(${up.length + 1}..${up.length + exon1.length},${e2s + 1}..${e2s + L2},${e3s + 1}..${e3s + exon3.length},${e4s + 1}..${e4s + exon4.length})`, '                     /gene="NR2F2"', '                     /transcript_id="NM_SYN3.1"', "ORIGIN",
  ...Array.from({ length: Math.ceil(genomic.length / 60) }, (_, i) => String(i * 60 + 1).padStart(9) + " " + genomic.slice(i * 60, i * 60 + 60).toLowerCase()), "//"].join("\n");
const sim = (n) => path.join(root, "fixtures", "sim", n);

test("command line: --genbank adds a functional knockout assessment to the table, the report and the provenance", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-gene-"));
  copyFileSync(sim("nr2f2_control.ab1"), path.join(dir, "NR2F2_WT.700.ab1")); copyFileSync(sim("nr2f2_ko_clone_hom_m1_edited.ab1"), path.join(dir, "NR2F2_A01_hom.700.ab1")); copyFileSync(sim("nr2f2_ko_clone_het_wt_m1_edited.ab1"), path.join(dir, "NR2F2_A02_het.700.ab1"));
  const gbPath = path.join(dir, "gene.gb"); writeFileSync(gbPath, gb); const out = path.join(dir, "out");
  assert.equal(await cli(["run", dir, "--design", path.join(root, "fixtures", "design", "nr2f2_2xha_ssodn.design.json"), "--workflow", "knockout", "--genbank", gbPath, "--regions", "kinase domain:20-60", "--epitope", "10-15", "--threads", "2", "--out", out]), 0);
  const results = JSON.parse(readFileSync(path.join(out, "results.json"), "utf8")); assert.equal(results.length, 2);
  results.forEach((r) => { assert.ok(r.functional, `${r.name}: functional assessment`); assert.ok(r.summary.functionalKoPct !== undefined); });
  const het = results.find((r) => /het/.test(r.name)); const hom = results.find((r) => /hom/.test(r.name));
  assert.ok(het.summary.functionalKoPct < hom.summary.functionalKoPct, "the heterozygous clone has a lower functional score than the homozygous one");
  assert.equal(het.functional.call.status, "not_knockout"); assert.equal(het.decision.decision, "reject");
  assert.equal(hom.functional.call.status, "predicted_knockout"); assert.equal(hom.decision.decision, "accept"); assert.ok(hom.summary.functionalKoPct > 90);
  const csv = readFileSync(path.join(out, "results.csv"), "utf8"); assert.match(csv.split("\n")[0], /functional_ko_pct/);
  const run = JSON.parse(readFileSync(path.join(out, "run.json"), "utf8")); assert.ok(run.geneRecord && /^[0-9a-f]{64}$/.test(run.geneRecord.sha256));
  assert.match(readFileSync(path.join(out, "report.html"), "utf8"), /Functional knockout assessment/);
  assert.ok(translate(exon1 + exon2 + exon3 + exon4).stopNt === (exon1 + exon2 + exon3 + exon4).length - 3, "the synthetic gene is a clean open reading frame");
});
