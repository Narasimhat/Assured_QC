import test from "node:test";
import assert from "node:assert/strict";
import { parseGenBank } from "../src/genbank.js";
import { reverseComplement } from "../src/seq.js";
import { buildCodingContext, assessAllele, assessFunction, translate, buildModels, describeContext } from "../src/coding.js";

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
// independent expectation: apply the edit to the coding sequence directly (edit inside one exon) and translate
const expectedPtc = (cdsOffset, delCount, ins) => { const m = cds.slice(0, cdsOffset) + ins + cds.slice(cdsOffset + delCount); return translate(m).stopNt; };

test("a GenBank record is read: sequence, joined CDS, complement strand", () => {
  const gb = parseGenBank(gbPlus); assert.equal(gb.sequence, genomic); const f = gb.features.find((x) => x.type === "CDS");
  assert.equal(f.location.strand, "+"); assert.equal(f.location.segments.length, 4); assert.deepEqual(f.location.segments[1], [start[1], ends[1]]); assert.equal(f.qualifiers.transcript_id[0], "NM_SYN1.1");
  const m = buildModels(gb, { gene: "SYN" }).models[0]; assert.equal(m.cdsSeq, cds); assert.equal(m.protLen, 169); assert.ok(m.startsWithAtg && m.endsWithStop);
  const gm = parseGenBank(gbMinus); const mm = buildModels(gm, {}).models[0]; assert.equal(mm.strand, "-"); assert.equal(mm.cdsSeq, cds, "the minus-strand record gives the same coding sequence");
});

test("a frameshift in an internal exon is predicted to trigger NMD, and the stop is where an independent translation puts it", () => {
  const r = buildCodingContext(gbPlus, { reference: win, gene: "SYN" }); assert.ok(r.ok, r.error); const ctx = r.context;
  const at = exon2InWin[0] + 60; const a = assessAllele(ctx, allele(win, at, at + 1));
  assert.equal(a.class, "lof_nmd", a.summary); assert.equal(a.lof, true);
  assert.equal(a.details.ptcNt, expectedPtc(100 + 60, 1, ""));
  assert.equal(a.details.deltaNt, -1); assert.ok(a.details.nmd.predicted);
  const ins = assessAllele(ctx, allele(win, at, at, "A")); assert.equal(ins.class, "lof_nmd"); assert.equal(ins.details.ptcNt, expectedPtc(160, 0, "A"));
});

test("the same edits are called the same on the opposite strand of the record and with the window reversed", () => {
  const at = exon2InWin[0] + 60; const plus = buildCodingContext(gbPlus, { reference: win }).context; const ref = assessAllele(plus, allele(win, at, at + 1));
  const minusGb = buildCodingContext(gbMinus, { reference: win }); assert.ok(minusGb.ok, minusGb.error);
  const a = assessAllele(minusGb.context, allele(win, at, at + 1)); assert.equal(a.class, ref.class); assert.equal(a.details.ptcNt, ref.details.ptcNt);
  const rwin = reverseComplement(win); const rc = buildCodingContext(gbPlus, { reference: rwin }); assert.ok(rc.ok, rc.error);
  const reversedAllele = reverseComplement(allele(win, at, at + 1)); const b = assessAllele(rc.context, reversedAllele);
  assert.equal(b.class, ref.class); assert.equal(b.details.ptcNt, ref.details.ptcNt);
  // a window with a few differences from the record (a control read with a SNP) is still located
  const snp = win.slice(0, 30) + (win[30] === "A" ? "C" : "A") + win.slice(31); assert.ok(buildCodingContext(gbPlus, { reference: snp }).ok);
});

test("in-frame deletions, the last coding exon, and edits outside the coding sequence", () => {
  const ctx = buildCodingContext(gbPlus, { reference: win }).context;
  const at = exon2InWin[0] + 59; const inframe = assessAllele(ctx, allele(win, at, at + 3)); // coding offset 159: one whole codon (residue 54)
  assert.equal(inframe.class, "inframe", inframe.summary); assert.equal(inframe.lof, null); assert.equal(inframe.details.diff.removed, 1); assert.equal(inframe.details.diff.first, 54);
  const withRegion = buildCodingContext(gbPlus, { reference: win, regions: [{ name: "kinase domain", start: 30, end: 60 }] }).context;
  assert.equal(assessAllele(withRegion, allele(win, at, at + 3)).class, "inframe_region_hit", "residue 54 is inside the marked region");
  const shifted = assessAllele(ctx, allele(win, at + 1, at + 4)); assert.equal(shifted.details.diff.removed + shifted.details.diff.inserted >= 1, true, "a deletion that straddles two codons replaces them");
  // last coding exon (100 nt): pick a one-base deletion whose new frame stops inside the exon; no NMD, and the protein keeps most of its length
  let off = 20; while (off < 90 && !(expectedPtc(410 + off, 1, "") >= 0)) off += 1;
  const w4 = genomic.slice(start[3] - 150, ends[3] + 50); const c4 = buildCodingContext(gbPlus, { reference: w4 }).context;
  const last = assessAllele(c4, allele(w4, 150 + off, 150 + off + 1)); assert.equal(last.details.nmd.predicted, false, last.summary); assert.equal(last.lof, null);
  assert.ok(last.details.retainedPct > 80, `retained ${last.details.retainedPct}`); assert.equal(last.details.ptcNt, expectedPtc(410 + off, 1, ""));
  const intronMid = introns[0].length >> 1; const w1i = genomic.slice(ends[0] - 30, start[1] + 30); const ci = buildCodingContext(gbPlus, { reference: w1i }).context;
  const out = assessAllele(ci, allele(w1i, 30 + intronMid, 30 + intronMid + 1)); assert.equal(out.class, "no_coding_effect"); assert.equal(out.lof, false);
  assert.equal(assessAllele(ctx, win).class, "wt");
});

test("a destroyed splice donor is modelled as exon skipping and intron retention", () => {
  const w = genomic.slice(ends[1] - 80, start[2] + 80); const ctx = buildCodingContext(gbPlus, { reference: w }).context;
  const donorAt = 80; const a = assessAllele(ctx, allele(w, donorAt, donorAt + 2)); // the GT after exon 2
  assert.ok(["lof_splice", "splice_uncertain"].includes(a.class), a.class); assert.ok(a.details.splice.skippedExons.includes("exon 2")); assert.ok(a.details.splice.skipping && a.details.splice.retention);
  // exon 2 is 150 nt: skipping it keeps the frame, so the call cannot be a firm loss of function
  assert.equal(a.details.splice.skipping.class, "inframe"); assert.equal(a.class, "splice_uncertain");
  // exon 3 is 160 nt: skipping it shifts the frame
  const w3 = genomic.slice(ends[2] - 80, start[3] + 80); const c3 = buildCodingContext(gbPlus, { reference: w3 }).context;
  const b = assessAllele(c3, allele(w3, 80, 82)); assert.notEqual(b.details.splice.skipping.class, "inframe");
  // the whole exon 2 with both flanking dinucleotides' neighbours deleted in one two-guide deletion: neighbours stay spliced, in frame
  const wb = genomic.slice(ends[0] - 60, start[2] + 60); const cb = buildCodingContext(gbPlus, { reference: wb }).context;
  const del = assessAllele(cb, allele(wb, 60 + 10, wb.length - 60 - 10 + 0)); // from inside intron 1 to inside intron 2, removing exon 2 and both dinucleotide sites
  assert.ok(del.class === "inframe" || del.class === "splice_uncertain" || del.class === "lof_splice", del.class);
});

test("start codon loss, and an isoform without the edited exon", () => {
  const w = genomic.slice(start[0] - 20, ends[0] + 20); const ctx = buildCodingContext(gbPlus, { reference: w }).context;
  const a = assessAllele(ctx, allele(w, 20, 23)); assert.ok(["lof_start_loss", "start_loss_uncertain"].includes(a.class), a.class);
  const two = buildCodingContext(gbTwo, { reference: win }); assert.equal(describeContext(two.context).isoforms.length, 2);
  const at = exon2InWin[0] + 60; const f = assessAllele(two.context, allele(win, at, at + 1)); assert.ok(f.details.isoforms.note, "one isoform has no exon 2"); assert.equal(f.details.isoforms.affected, 1);
});

test("the clone call: knockout, heterozygous, in frame", () => {
  const ctx = buildCodingContext(gbPlus, { reference: win }).context; const at = exon2InWin[0] + 60;
  const fs1 = { id: "a", kind: "indel", label: "-1", size: -1, sequence: allele(win, at, at + 1) }; const fs2 = { id: "b", kind: "indel", label: "+1", size: 1, sequence: allele(win, at, at, "G") };
  const inf = { id: "c", kind: "indel", label: "-3", size: -3, sequence: allele(win, at + 9, at + 12) }; const wt = { id: "wt", kind: "wt", label: "Wild type", size: 0, sequence: win };
  const hom = assessFunction(ctx, [{ ...fs1, fraction: 0.97 }, { ...wt, fraction: 0.03 }]); assert.equal(hom.call.status, "predicted_knockout"); assert.ok(hom.score.lofPct >= 96);
  const compound = assessFunction(ctx, [{ ...fs1, fraction: 0.5 }, { ...fs2, fraction: 0.5 }]); assert.equal(compound.call.status, "predicted_knockout"); assert.equal(compound.score.lofPct, 100);
  const het = assessFunction(ctx, [{ ...fs1, fraction: 0.5 }, { ...wt, fraction: 0.5 }]); assert.equal(het.call.status, "not_knockout"); assert.equal(het.score.lofPct, 50);
  const mixed = assessFunction(ctx, [{ ...fs1, fraction: 0.5 }, { ...inf, fraction: 0.5 }]); assert.equal(mixed.call.status, "uncertain"); assert.equal(mixed.score.uncertainPct, 50);
  assert.ok(compound.confirm.some((s) => /Western/.test(s)) && compound.confirm.some((s) => /RT-qPCR/.test(s)));
});

test("a deletion that can be written two ways is not called a splice-site loss when the donor GT is re-formed at the junction", () => {
  // find an n so that the last base of exon 3 minus n equals the donor G: deleting [exon end - n + 1, donor G] is the same allele as deleting [exon end - n, exon end)
  const w = genomic.slice(ends[2] - 80, start[3] + 80); const ctx = buildCodingContext(gbPlus, { reference: w }).context; const donor = 80;
  let n = 0; for (let k = 4; k < 40; k += 1) if (w[donor - k] === w[donor] && w[donor - k + 1] !== w[donor + 1]) { n = k; break; }
  assert.ok(n > 0, "found a repeat of the donor base inside the exon");
  const right = allele(w, donor - n + 1, donor + 1); const left = allele(w, donor - n, donor); assert.equal(right, left, "the two placements give one allele");
  const a = assessAllele(ctx, right);
  assert.ok(!["lof_splice", "splice_uncertain"].includes(a.class), `read as an exon shortened by ${n} nt, not as a lost donor: ${a.class}`); assert.equal(a.details.spliceReformed, true); assert.match(a.summary, /re-formed|consensus/);
  // a deletion that cannot be shifted off the donor still loses it
  let m = 0; for (let k = 2; k < 30; k += 1) if (w[donor - k] !== w[donor] && w[donor - k + 1] !== w[donor + 1] && w[donor - k] !== w[donor + 1]) { m = k; break; }
  const hit = assessAllele(ctx, allele(w, donor - m + 1, donor + 1)); assert.ok(["lof_splice", "splice_uncertain"].includes(hit.class), hit.class);
});

test("a record that also carries a neighbouring gene: only the gene the reference window lies in is modelled", () => {
  const w = genomic.slice(start[1] + 20, start[1] + 140);
  // a longer neighbouring gene on each strand, appended after the design's gene
  const neighbour = "ATG" + Array.from({ length: 400 }, () => CODONS[Math.floor(rnd() * CODONS.length)]).join("") + "TAA";
  const pad = randomDna(200); const seq = genomic + pad + neighbour + randomDna(200); const s0 = genomic.length + pad.length;
  const rec = (strand) => genBank(seq, [{ type: "CDS", loc: joinText, q: { gene: "SYN", transcript_id: "NM_SYN1.1" } }, { type: "CDS", loc: strand === "+" ? `${s0 + 1}..${s0 + neighbour.length}` : `complement(${s0 + 1}..${s0 + neighbour.length})`, q: { gene: "OTHER", transcript_id: "NM_OTHER.1" } }]);
  for (const strand of ["+", "-"]) {
    const built = buildCodingContext(rec(strand), { reference: w }); assert.ok(built.ok, built.error);
    assert.equal(built.context.models.length, 1, `the ${strand}-strand neighbour is left out`); assert.equal(built.context.models[0].gene, "SYN"); assert.ok(built.context.warnings.some((x) => /other genes/.test(x)));
    assert.equal(built.context.models[0].protLen, 169);
  }
  // a window that lies in no gene keeps every coding sequence
  const outside = buildCodingContext(rec("+"), { reference: genomic.slice(0, 120) }); assert.ok(outside.ok); assert.ok(outside.context.models.length >= 1);
});
