import test from "node:test";
import assert from "node:assert/strict";
import { mkdtempSync, mkdirSync, copyFileSync, writeFileSync, readFileSync, existsSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { deflateRawSync } from "node:zlib";
import { readZip, readIceExport, parseIceSummary, parseIndels, looksLikeIceSummary, matchIceRows, compareWithExport } from "../src/iceImport.js";
import { iceSectionHtml, iceCsv } from "../src/iceCompare.js";
import { main as cli } from "../cli/assured-qc.mjs";

const here = path.dirname(fileURLToPath(import.meta.url)); const root = path.join(here, "..");

// the summary of an ICE web export (layout of ICE 4.2); sample names, guides and numbers here are made up
const G1 = "ACGTACGTACGTACGTACGT", G2 = "TTGACCAGTTGACCAGTTGA";
const HEADER = "Label,ICE,KO-Score,KI-Score,ICE d,R Squared,Mean Discord Before,Mean Discord After,Guide Sequences,Control Sample Quality Score,Edit Sample Quality Score,Notes,Indels";
const row = (label, ice, ko, ki, r2, indels, { guides = `${G1},${G2}`, notes = "" } = {}) => `${label},${ice},${ko},${ki},${ice},${r2},0.13,0.72,"${guides}",58,60,${notes},"${indels}"`;
const SUMMARY = [HEADER,
  row("p1_gene_KO", 100, 100, "None", 1.0, "{'-78': 100.0}"),
  row("p2_gene_KO", 100, 100, "None", 0.96, "{'-77': 49.0, '-79': 51.0}"),
  row("p3_gene_KO", 58, 58, "None", 0.63, "{'-78': 41.0, '+5': 5.0, '+4': 3.0}", { notes: "low fit" }),
  row("p4_gene_KO", 100, 100, "None", 1.0, "{'-79': 100.0}", { guides: G1 }),
  row("p10_gene_KO", 100, 100, "None", 1.0, "{'-79': 100.0}"),
  row("other_plate_7", 100, 100, "None", 1.0, "{'-79': 100.0}")].join("\n") + "\n";

// a zip writer for the tests: stored or deflated entries, central directory, end record
const crcTable = Array.from({ length: 256 }, (_, n) => { let c = n; for (let k = 0; k < 8; k += 1) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; return c >>> 0; });
const crc32 = (b) => { let c = 0xffffffff; for (const x of b) c = crcTable[(c ^ x) & 0xff] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; };
function makeZip(entries) { // entries: [name, string, "stored" | "deflate"]
  const chunks = [], central = []; let offset = 0;
  const w16 = (n) => Buffer.from([n & 255, (n >> 8) & 255]); const w32 = (n) => Buffer.from([n & 255, (n >> 8) & 255, (n >> 16) & 255, (n >>> 24) & 255]);
  for (const [name, text, how = "deflate"] of entries) {
    const raw = Buffer.from(text); const data = how === "stored" || !raw.length ? raw : deflateRawSync(raw); const method = how === "stored" || !raw.length ? 0 : 8; const nm = Buffer.from(name);
    const local = Buffer.concat([w32(0x04034b50), w16(20), w16(0), w16(method), w16(0), w16(0), w32(crc32(raw)), w32(data.length), w32(raw.length), w16(nm.length), w16(0), nm, data]);
    central.push(Buffer.concat([w32(0x02014b50), w16(20), w16(20), w16(0), w16(method), w16(0), w16(0), w32(crc32(raw)), w32(data.length), w32(raw.length), w16(nm.length), w16(0), w16(0), w16(0), w16(0), w32(0), w32(offset), nm]));
    chunks.push(local); offset += local.length;
  }
  const dir = Buffer.concat(central);
  return new Uint8Array(Buffer.concat([...chunks, dir, w32(0x06054b50), w16(0), w16(0), w16(entries.length), w16(entries.length), w32(dir.length), w32(offset), w16(0)]));
}

test("the Indels cell of the summary is read as sizes and percentages; names that are not sizes are left out", () => {
  assert.deepEqual(parseIndels("{'-77': 49.0, '-79': 51.0}"), [{ size: -77, pct: 49 }, { size: -79, pct: 51 }]);
  assert.deepEqual(parseIndels("{'+5': 5.0, '0': 2.5, 'HDR': 80.0}"), [{ size: 5, pct: 5 }, { size: 0, pct: 2.5 }]);
  assert.deepEqual(parseIndels(""), []);
});

test("the summary is read: scores, None as no value, guides, alleles; a table that is not an ICE summary is refused", () => {
  const { rows, problems } = parseIceSummary(SUMMARY); assert.deepEqual(problems, []); assert.equal(rows.length, 6);
  assert.deepEqual(rows[1], { label: "p2_gene_KO", ice: 100, ko: 100, ki: null, rsq: 0.96, guides: [G1, G2], notes: "", indels: [{ size: -77, pct: 49 }, { size: -79, pct: 51 }] });
  assert.equal(rows[2].notes, "low fit"); assert.deepEqual(rows[3].guides, [G1]);
  assert.match(parseIceSummary("a,b\n1,2\n").problems[0], /no Label column/); assert.match(parseIceSummary("Label,Name\nx,y\n").problems[0], /no ICE, KO-Score or KI-Score/);
  assert.equal(looksLikeIceSummary(SUMMARY), true); assert.equal(looksLikeIceSummary("sample,control,edited\na,b,c\n"), false);
});

test("labels are matched to samples by whole words: cl4 is not in cl40, the longest name wins, an equal fit is not assigned, what is left is the donor", () => {
  const names = ["cl4", "cl40", "A7", "A7 b", "x1", "X1"];
  const m = matchIceRows(names, [{ label: "cl4_KO_run2" }, { label: "cl40_KO" }, { label: "A7_b_ssODN1.ab1" }, { label: "A7-ssODN2" }, { label: "x1_gene" }, { label: "nothing" }]);
  assert.deepEqual([...m.bySample.keys()].sort(), [0, 1, 2, 3]);
  assert.equal(m.bySample.get(0)[0].row.label, "cl4_KO_run2"); assert.equal(m.bySample.get(0).length, 1); assert.equal(m.bySample.get(1)[0].rest, "KO");
  assert.equal(m.bySample.get(3)[0].rest, "ssODN1", "the longer name 'A7 b' wins over 'A7'"); assert.equal(m.bySample.get(2)[0].rest, "ssODN2");
  assert.deepEqual(m.ambiguous.map((a) => a.row.label), ["x1_gene"]); assert.deepEqual(m.unmatched.map((r) => r.label), ["nothing"]);
});

test("a zip is read: stored and deflated entries, folders skipped, the top-level summary chosen over the one in a sample folder; damaged and foreign files give a message", async () => {
  const zip = makeZip([["results/", ""], ["summary_abc.csv", SUMMARY], ["results/p1_gene_KO/summary_old.csv", "Label,KO-Score,R Squared\nzz,1,1\n", "stored"], ["results/p1_gene_KO/all.txt", "x".repeat(5000)], ["Understanding your ICE Export.pdf", "%PDF", "stored"]]);
  const entries = await readZip(zip, (n) => /summary/.test(n)); assert.deepEqual([...entries.keys()].sort(), ["results/p1_gene_KO/summary_old.csv", "summary_abc.csv"]);
  assert.equal(new TextDecoder().decode(entries.get("summary_abc.csv")), SUMMARY);
  const read = await readIceExport([{ name: "ICE_Results_x.zip", bytes: zip }]); assert.deepEqual(read.problems, []); assert.equal(read.rows.length, 6); assert.deepEqual(read.sources, ["ICE_Results_x.zip"]);
  const csv = await readIceExport([{ name: "summary_abc.csv", bytes: new TextEncoder().encode(SUMMARY) }]); assert.equal(csv.rows.length, 6);
  const none = await readIceExport([{ name: "empty.zip", bytes: makeZip([["notes.txt", "hello"]]) }]); assert.match(none.problems[0], /no summary \.csv in the zip/);
  const damaged = await readIceExport([{ name: "bad.zip", bytes: zip.slice(0, zip.length - 30) }]); assert.equal(damaged.rows.length, 0); assert.ok(damaged.problems.length);
  const wrong = await readIceExport([{ name: "sheet.csv", text: "sample,control,edited\na,b,c\n" }]); assert.match(wrong.problems[0], /no ICE, KO-Score or KI-Score column/);
});

// results as a run makes them: knockout clones with two deletions of 77 and 79 bases
const ko = (name, sizes = [[-77, 0.5], [-79, 0.5]]) => ({ ok: true, name, workflow: "deletion", sampleType: "clone", summary: { koScorePct: 100, wtPct: 0 }, intervals: { koScorePct: [97, 100] }, contributions: sizes.map(([size, fraction]) => ({ kind: "deletion_between_cuts", size, fraction })) });
const task = (name, workflow = "deletion", guides = [G1, G2], extra = {}) => ({ name, workflow, spec: { guides: guides.map((spacer) => ({ spacer })) }, controls: ["wt.ab1"], edited: [`${name}.ab1`], ...extra });

test("a knockout plate against an ICE export: agree, a poor fit that is not a difference, other guides flagged, samples without a result, labels without a sample", () => {
  const tasks = ["p1", "p2", "p3", "p4", "p5"].map((n) => task(n)); const results = [ko("p1", [[-78, 1]]), ko("p2"), ko("p3", [[-78, 0.5], [-77, 0.5]]), ko("p4", [[-79, 1]]), ko("p5")];
  const read = parseIceSummary(SUMMARY); const cmp = compareWithExport({ tasks, results, rows: read.rows });
  const by = Object.fromEntries(cmp.rows.map((r) => [r.sample, r]));
  assert.equal(by.p1.status, "agree"); assert.equal(by.p1.ice, 100); assert.equal(by.p2.status, "agree"); assert.deepEqual(by.p2.iceAlleles.map((a) => a.size).sort((a, b) => a - b), [-79, -77]);
  assert.equal(by.p3.status, "ice_fit_poor"); assert.equal(by.p3.iceR2, 0.63); assert.match(by.p3.note, /Assured QC also reads -77 \(50%\), which ICE does not/); assert.match(by.p3.note, /ICE note: low fit/);
  assert.equal(by.p4.status, "agree"); assert.match(by.p4.note, /ICE was run with 1 guide/);
  assert.equal(by.p5.status, "not_run"); assert.match(by.p5.note, /no sample of the ICE export has this sample's name/);
  assert.equal(cmp.unmatched.map((r) => r.label).join(), "p10_gene_KO,other_plate_7", "p10 is not p1");
  assert.match(cmp.problems[0], /2 labels of the ICE export match none of the run's samples: p10_gene_KO, other_plate_7/);
});

test("a point edit with one ICE run per donor: the donor comes from the rest of the label, the KI-Score is used and the best donor is shown", () => {
  const rows = parseIceSummary([HEADER, row("a1_ssODN1", 70, 20, 62, 0.97, "{'-3': 8.0}", { guides: G1 }), row("a1_ssODN2", 80, 30, 81, 0.98, "{'-3': 8.0}", { guides: G1 })].join("\n")).rows;
  const edit = { ok: true, name: "a1", workflow: "snp", sampleType: "clone", summary: { intendedEditPct: 85, wtPct: 5 }, intervals: { intendedEditPct: [78, 100] }, contributions: [] };
  const cmp = compareWithExport({ tasks: [task("a1", "snp", [G1])], results: [edit], rows });
  assert.equal(cmp.rows.length, 1); const r = cmp.rows[0]; assert.equal(r.status, "agree"); assert.equal(r.ice, 81); assert.equal(r.donor, "ssODN2"); assert.match(r.note, /Best of 2 donors run one at a time \(ssODN1: 62, ssODN2: 81\)/);
});

test("the report section says that the results come from the web tool and that the ICE version matters; the CSV carries the same rows", () => {
  const cmp = compareWithExport({ tasks: [task("p1")], results: [ko("p1", [[-78, 1]])], rows: parseIceSummary(SUMMARY).rows });
  const html = iceSectionHtml(cmp.rows, { imported: true, source: "ICE_Results_x.zip" });
  assert.match(html, /Results exported from the ICE web tool \(ICE_Results_x\.zip\)/); assert.match(html, /ICE version matters/); assert.doesNotMatch(html, /non-profit/);
  assert.match(iceCsv(cmp.rows), /^sample,workflow,assured_qc_score/); assert.match(iceSectionHtml(cmp.rows, { version: "1.2.0", python: "3.11" }), /ICE 1\.2\.0 run separately/);
});

test("command line: --ice-results reads the export, writes the comparison and the provenance, and refuses a missing file or a combination with --ice", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-ice-import-")); const sim = (n) => path.join(root, "fixtures", "sim", n);
  const files = { "apoe_WT.700.ab1": "apoe_control.ab1", "apoe_A01_hom.700.ab1": "apoe_snp_clone_hom_edited.ab1", "apoe_A02_het.700.ab1": "apoe_snp_clone_het_edited.ab1" };
  Object.entries(files).forEach(([to, from]) => copyFileSync(sim(from), path.join(dir, to)));
  const design = path.join(root, "fixtures", "design", "apoe_r176c_snp.design.json"); const out = path.join(dir, "out");
  const summary = [HEADER, row("apoe_A01_hom.700_ssODN1", 95, 5, 94, 0.99, "{}", { guides: "" }), row("apoe_A02_het.700_ssODN1", 60, 5, 20, 0.95, "{}", { guides: "" }), row("not_in_run", 100, 100, "None", 1, "{'-5': 100.0}", { guides: "" })].join("\n") + "\n";
  const zipPath = path.join(dir, "ICE_Results_test.zip"); writeFileSync(zipPath, makeZip([["summary_test.csv", summary], ["results/x/all.txt", "x"]]));
  assert.equal(await cli(["run", dir, "--design", design, "--threads", "2", "--out", out, "--ice-results", zipPath]), 0);
  const csv = readFileSync(path.join(out, "ice_comparison.csv"), "utf8").trim().split("\n"); assert.equal(csv.length, 3, "header and the two samples");
  const run = JSON.parse(readFileSync(path.join(out, "run.json"), "utf8")); assert.equal(run.ice.imported, true); assert.equal(run.ice.files[0].file, "ICE_Results_test.zip"); assert.match(run.ice.files[0].sha256, /^[0-9a-f]{64}$/); assert.deepEqual(run.ice.labelsNotMatched, ["not_in_run"]);
  assert.match(readFileSync(path.join(out, "report.html"), "utf8"), /Results exported from the ICE web tool/); assert.ok(existsSync(path.join(out, "ice_results_used.json")));
  const used = JSON.parse(readFileSync(path.join(out, "ice_results_used.json"), "utf8")); assert.equal(used.rows.length, 3);
  assert.equal(await cli(["run", dir, "--design", design, "--out", path.join(dir, "o2"), "--ice-results", path.join(dir, "missing.zip")]), 1);
  assert.equal(await cli(["run", dir, "--design", design, "--out", path.join(dir, "o3"), "--ice-results", zipPath, "--ice", dir]), 1);
  // the same results from the unzipped summary and from a folder that holds it
  const folder = path.join(dir, "export"); mkdirSync(folder); writeFileSync(path.join(folder, "summary_test.csv"), summary);
  assert.equal(await cli(["run", dir, "--design", design, "--threads", "2", "--out", path.join(dir, "o4"), "--ice-results", folder]), 0);
  assert.equal(readFileSync(path.join(dir, "o4", "ice_comparison.csv"), "utf8"), readFileSync(path.join(out, "ice_comparison.csv"), "utf8"));
});
