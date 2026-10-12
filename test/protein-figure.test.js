import test from "node:test";
import assert from "node:assert/strict";
import { reverseComplement } from "../src/seq.js";
import { buildCodingContext, assessFunction, translate } from "../src/coding.js";
import { predictedProducts, epitopeWindows, drawProteinProducts, cloneProteinFigure, runProteinFigure, joinAlleles } from "../src/proteinFigure.js";

// A synthetic four-exon gene, + strand. Coding lengths 100 / 150 / 160 / 100 (510 nt, 169 residues + stop); the exon junctions fall in different phases.
let seed = 7; const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
const BASES = "ACGT";
const CODONS = []; for (const a of BASES) for (const b of BASES) for (const c of BASES) if (!["TAA", "TAG", "TGA"].includes(a + b + c) && a + b + c !== "ATG") CODONS.push(a + b + c);
const randomDna = (n) => Array.from({ length: n }, () => BASES[Math.floor(rnd() * 4)]).join("");
const cds = "ATG" + Array.from({ length: 168 }, () => CODONS[Math.floor(rnd() * CODONS.length)]).join("") + "TAA";
const exonCds = [cds.slice(0, 100), cds.slice(100, 250), cds.slice(250, 410), cds.slice(410)];
const utr5 = randomDna(40); const up = randomDna(300); const down = randomDna(300);
const introns = [0, 1, 2].map(() => "GT" + randomDna(296) + "AG");
const genomic = up + utr5 + exonCds[0] + introns[0] + exonCds[1] + introns[1] + exonCds[2] + introns[2] + exonCds[3] + down;
const start = [up.length + utr5.length]; // genomic start of each exon's coding part (0-based)
start.push(start[0] + 100 + 300, 0, 0); start[2] = start[1] + 150 + 300; start[3] = start[2] + 160 + 300;
const ends = [start[0] + 100, start[1] + 150, start[2] + 160, start[3] + 100];
assert.equal(genomic.slice(start[1], ends[1]), exonCds[1]);

function genBank(sequence, features) {
  const lines = ["LOCUS       SYN_GENE   " + sequence.length + " bp    DNA     linear   SYN", "DEFINITION  synthetic four-exon gene.", "ACCESSION   NG_SYN1", "FEATURES             Location/Qualifiers"];
  features.forEach((f) => { lines.push(`     ${f.type.padEnd(16)}${f.loc}`); Object.entries(f.q || {}).forEach(([k, v]) => lines.push(`                     /${k}="${v}"`)); });
  lines.push("ORIGIN"); for (let i = 0; i < sequence.length; i += 60) lines.push(String(i + 1).padStart(9) + " " + sequence.slice(i, i + 60).toLowerCase().replace(/(.{10})/g, "$1 "));
  lines.push("//"); return lines.join("\n");
}
const segs = start.map((s, i) => [s + 1, ends[i]]); // 1-based inclusive
const joinText = `join(${segs.map(([a, b]) => `${a}..${b}`).join(",")})`;
const isoJoin = `join(${[segs[0], segs[2], segs[3]].map(([a, b]) => `${a}..${b}`).join(",")})`;
const gbPlus = genBank(genomic, [{ type: "gene", loc: `${segs[0][0] - 40}..${segs[3][1]}`, q: { gene: "SYN" } }, { type: "CDS", loc: joinText, q: { gene: "SYN", transcript_id: "NM_SYN1.1", protein_id: "NP_SYN1.1", product: "synthetic protein" } }]);
const L = genomic.length; const flipLoc = (join) => `complement(join(${join.replace(/join\(|\)/g, "").split(",").map((p) => { const [a, b] = p.split("..").map(Number); return `${L - b + 1}..${L - a + 1}`; }).reverse().join(",")}))`;
const gbMinus = genBank(reverseComplement(genomic), [{ type: "CDS", loc: flipLoc(joinText), q: { gene: "SYN", transcript_id: "NM_SYN1.1", protein_id: "NP_SYN1.1" } }]);
const gbTwo = genBank(genomic, [{ type: "CDS", loc: joinText, q: { gene: "SYN", transcript_id: "NM_SYN1.1" } }, { type: "CDS", loc: isoJoin, q: { gene: "SYN", transcript_id: "NM_SYN2.1" } }]);

// reference windows (what the QC app holds): from inside exon 1 to inside exon 4
const w0 = start[1] - 120; const w1 = ends[2] + 120; const win = genomic.slice(w0, w1);
const cdsOf = (g) => [[start[0], ends[0]], [start[1], ends[1]], [start[2], ends[2]], [start[3], ends[3]]].map(([a, b]) => g.slice(a, b)).join("");
// helper: an allele window made by replacing window[s, e) with ins
const allele = (W, s, e, ins = "") => W.slice(0, s) + ins + W.slice(e);
const exon2InWin = [start[1] - w0, ends[1] - w0]; const exon3InWin = [start[2] - w0, ends[2] - w0];

const fnOf = (ctx, alleles) => assessFunction(ctx, alleles.map((a, i) => ({ id: `a${i}`, kind: "indel", label: a.label, size: a.size ?? -1, fraction: a.fraction, sequence: a.sequence })));
const result = (functional, name = "clone1") => ({ ok: true, name, kind: "trace-pair", workflow: "knockout", functional, summary: {}, warnings: [] });

test("a frameshift gives the residues kept before the change, a frameshifted tail to the stop, and the epitope windows", () => {
  const ctx = buildCodingContext(gbPlus, { reference: win, gene: "SYN" }).context; const at = exon2InWin[0] + 60;
  const fn = fnOf(ctx, [{ label: "-1 bp", sequence: allele(win, at, at + 1), fraction: 1 }]);
  const data = predictedProducts(fn); assert.equal(data.proteinLength, 169); assert.equal(data.exonEnds.length, 4);
  assert.equal(data.products.length, 1); const p = data.products[0];
  assert.equal(p.kind, "truncated"); assert.equal(p.total, fn.rows[0].details.ptcAa); assert.equal(p.kept, fn.rows[0].details.firstChangedAa - 1); assert.ok(p.total >= p.kept && p.kept > 40);
  const w = epitopeWindows(data.products, data.proteinLength); assert.deepEqual(w.all, [[1, p.kept]]); assert.deepEqual(w.only, [[p.kept + 1, 169]]);
  const svg = cloneProteinFigure(result(fn)); assert.match(svg, /^<svg /); assert.match(svg, new RegExp(`${p.total} aa, about ${(Math.round(p.total * 0.11 * 10) / 10)} kDa`)); assert.match(svg, new RegExp(`residues 1-${p.kept} sees every product`));
  assert.match(svg, /Full-length protein: 169 aa/); assert.ok(!/<script|javascript:/i.test(svg)); assert.equal((svg.match(/<svg /g) || []).length, 1); assert.ok(svg.trimEnd().endsWith("</svg>"));
});

test("a splice-donor hit gives exon skipping (the residues before the skipped exon are kept) and intron retention", () => {
  // the donor after exon 3 (160 nt, so skipping it shifts the frame): 100 + 150 nt before it = 83 complete codons
  const w = genomic.slice(ends[2] - 80, start[3] + 80); const ctx = buildCodingContext(gbPlus, { reference: w, gene: "SYN" }).context;
  const fn = fnOf(ctx, [{ label: "donor hit", size: -2, sequence: allele(w, 80, 82), fraction: 1 }]); assert.ok(["lof_splice", "splice_uncertain"].includes(fn.rows[0].class), fn.rows[0].class);
  const data = predictedProducts(fn); const skip = data.products.find((p) => /exon skipping/.test(p.outcome)); const keep = data.products.find((p) => p.outcome === "intron retained");
  assert.ok(skip && keep, JSON.stringify(data.products.map((p) => p.outcome)));
  assert.match(skip.outcome, /exon 3/); assert.equal(skip.kind, "truncated"); assert.equal(skip.kept, 83, "exons 1 and 2 are 250 nt: 83 complete codons before the skipped exon");
  assert.equal(skip.total, fn.rows[0].details.splice.skipping.details.ptcAa); assert.equal(keep.kind, "truncated"); assert.equal(keep.total, fn.rows[0].details.splice.retention.details.ptcAa);
  const w2 = epitopeWindows(data.products, 169); assert.ok(w2.all.length && w2.all[0][1] <= Math.min(skip.kept, keep.kept));
  const svg = drawProteinProducts(data, { title: "t" }); assert.match(svg, /exon skipping \(exon 3\)/); assert.match(svg, /intron retained/);
  // an exon of 150 nt (50 codons) skipped in frame is an in-frame loss of those residues, drawn as removed
  const w1 = genomic.slice(ends[1] - 80, start[2] + 80); const c1 = buildCodingContext(gbPlus, { reference: w1, gene: "SYN" }).context;
  const inframe = predictedProducts(fnOf(c1, [{ label: "donor hit", size: -2, sequence: allele(w1, 80, 82), fraction: 1 }])).products.find((p) => /exon skipping/.test(p.outcome));
  assert.equal(inframe.kind, "inframe"); assert.ok(inframe.removed === 50 || inframe.removed === 51, "50 codons, plus the codon that joins exon 1 to exon 3 when its phase changes"); assert.equal(inframe.kept, 33);
});

test("an in-frame deletion is drawn with the removed residues; a start-codon hit is listed as text; synonymous and wild-type alleles draw nothing", () => {
  const ctx = buildCodingContext(gbPlus, { reference: win, gene: "SYN" }).context; const at = exon2InWin[0] + 59;
  const fn = fnOf(ctx, [{ label: "-3 bp", size: -3, sequence: allele(win, at, at + 3), fraction: 0.5 }, { label: "WT", sequence: win, fraction: 0.5 }]);
  const data = predictedProducts(fn); assert.equal(data.products.length, 1); const p = data.products[0];
  assert.equal(p.kind, "inframe"); assert.equal(p.removed, 1); assert.equal(p.kept, 53); assert.equal(p.resume, 55); assert.equal(p.total, 168);
  const svg = drawProteinProducts(data); assert.match(svg, /1 residue removed at residue 54/); assert.match(svg, /stroke-dasharray="2 2"/);
  const w = genomic.slice(start[0] - 20, ends[0] + 20); const c0 = buildCodingContext(gbPlus, { reference: w, gene: "SYN" }).context;
  const start_ = predictedProducts(fnOf(c0, [{ label: "start hit", size: -3, sequence: allele(w, 20, 23), fraction: 1 }])); assert.ok(start_.products.some((x) => x.kind === "text"), "an allele without a drawable product is listed as text");
  assert.equal(epitopeWindows(start_.products, 169), null); assert.match(drawProteinProducts(start_), /start hit/);
  const wtOnly = fnOf(ctx, [{ label: "WT", sequence: win, fraction: 1 }]); assert.equal(predictedProducts(wtOnly).products.length, 0); assert.equal(cloneProteinFigure(result(wtOnly)), "");
});

test("marked regions and the user's epitope are drawn; products with the same outcome from different alleles are one row", () => {
  const ctx = buildCodingContext(gbPlus, { reference: win, gene: "SYN", regions: [{ name: "kinase domain", start: 30, end: 120 }], epitope: { start: 10, end: 25 } }).context; const at = exon2InWin[0] + 60;
  const fn = fnOf(ctx, [{ label: "-1 bp (a)", sequence: allele(win, at, at + 1), fraction: 0.5 }, { label: "-1 bp (b)", sequence: allele(win, at + 1, at + 2, ""), fraction: 0.5 }]);
  const data = predictedProducts(fn); assert.equal(data.regions.length, 1); assert.deepEqual(data.epitope, { start: 10, end: 25 });
  const svg = drawProteinProducts(data); assert.match(svg, /kinase domain/); assert.match(svg, /Your antibody epitope \(residues 10-25, dashed box\) is kept by \d of \d products/);
  assert.equal(joinAlleles(["Deletion between cuts (79 bp)", "Deletion between cuts (84 bp)"]), "Deletion between cuts (79 bp, 84 bp)"); assert.equal(joinAlleles(["a", "b"]), "a; b"); assert.equal(joinAlleles(["only (1 bp)"]), "only (1 bp)");
  const same = predictedProducts({ context: { proteinLength: 100, exonLengths: [150, 150] }, rows: [
    { label: "Del (3 bp)", class: "lof_nmd", lof: true, fraction: 0.5, details: { ptcAa: 40, firstChangedAa: 31, nmd: { predicted: true } } }, { label: "Del (4 bp)", class: "lof_nmd", lof: true, fraction: 0.5, details: { ptcAa: 40, firstChangedAa: 31, nmd: { predicted: true } } }] });
  assert.equal(same.products.length, 1); assert.equal(same.products[0].allele, "Del (3 bp, 4 bp)"); assert.ok(Math.abs(same.products[0].share - 1) < 1e-9);
});

test("text is escaped, nothing is drawn without a functional assessment, and a run figure merges the clones that share a protein", () => {
  const evil = predictedProducts({ context: { proteinLength: 100, exonLengths: [150, 150] }, rows: [{ label: "<b>x</b> & \"y\"", class: "lof_nmd", lof: true, fraction: 1, details: { ptcAa: 40, firstChangedAa: 31, nmd: { predicted: true } } }] });
  const svg = drawProteinProducts(evil, { title: "<img src=x onerror=1>" }); assert.ok(!/<img|<b>/.test(svg)); assert.match(svg, /&lt;b&gt;x&lt;\/b&gt; &amp; &quot;y&quot;/);
  assert.equal(cloneProteinFigure({ ok: true, name: "x" }), ""); assert.equal(cloneProteinFigure({ ok: false, error: "e" }), ""); assert.equal(runProteinFigure([{ ok: true, name: "x" }]), ""); assert.equal(runProteinFigure([]), "");
  const ctx = buildCodingContext(gbPlus, { reference: win, gene: "SYN" }).context; const at = exon2InWin[0] + 60;
  const f1 = fnOf(ctx, [{ label: "-1 bp", sequence: allele(win, at, at + 1), fraction: 1 }]); const f2 = fnOf(ctx, [{ label: "-1 bp", sequence: allele(win, at, at + 1), fraction: 1 }, { label: "-2 bp", size: -2, sequence: allele(win, at, at + 2), fraction: 0 }]);
  const run = runProteinFigure([result(f1, "c1"), result(f2, "c2")]); assert.match(run, /alleles in this run \(2 clones\)/); assert.equal((run.match(/aa, about/g) || []).length, 2, "the full-length protein and one product");
});
