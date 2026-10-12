// A synthetic four-exon gene built around the NR2F2 reference window of the design fixture: exon 2 carries the cut, so a 1-nt deletion there puts a stop
// well upstream of the last exon junction (NMD expected). Shared by the command-line and design-file tests.
import { readFileSync } from "node:fs";
import path from "node:path";
import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";
import { translate } from "../../src/coding.js";

const here = path.dirname(fileURLToPath(import.meta.url)); const root = path.join(here, "..", "..");
export const spec = JSON.parse(readFileSync(path.join(root, "fixtures", "design", "nr2f2_2xha_ssodn.design.json"), "utf8")); const ref = spec.reference.sequence; const cut = spec.guides[0].cut;
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
export const gb = ["LOCUS       SYN3   " + genomic.length + " bp", "ACCESSION   NG_SYN3", "FEATURES             Location/Qualifiers", `     CDS             join(${up.length + 1}..${up.length + exon1.length},${e2s + 1}..${e2s + L2},${e3s + 1}..${e3s + exon3.length},${e4s + 1}..${e4s + exon4.length})`, '                     /gene="NR2F2"', '                     /transcript_id="NM_SYN3.1"', "ORIGIN",
  ...Array.from({ length: Math.ceil(genomic.length / 60) }, (_, i) => String(i * 60 + 1).padStart(9) + " " + genomic.slice(i * 60, i * 60 + 60).toLowerCase()), "//"].join("\n");

export const geneParts = { exon1, exon2, exon3, exon4, up, i1, i2, i3, a2, L2 };
export const windowStartInGene = up.length + exon1.length + i1.length; // where the design reference starts inside the synthetic gene
export const geneGenomic = genomic;
export const codingSegments = [[up.length, up.length + exon1.length], [e2s, e2s + L2], [e3s, e3s + exon3.length], [e4s, e4s + exon4.length]];
export const designRoot = root;
assert.ok(translate(exon1 + exon2 + exon3 + exon4).stopNt === (exon1 + exon2 + exon3 + exon4).length - 3, "the synthetic gene is a clean open reading frame");
