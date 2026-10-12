import test from "node:test";
import assert from "node:assert/strict";
import { mkdtempSync, copyFileSync, readFileSync, existsSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { main as cli } from "../cli/assured-qc.mjs";
import { decide, tally, wellOf } from "../src/decision.js";
import { buildPlateReportHtml, resequenceListCsv } from "../src/plateReport.js";

const here = path.dirname(fileURLToPath(import.meta.url)); const root = path.join(here, "..");
const sim = (n) => path.join(root, "fixtures", "sim", n);

test("command line: a folder, a design file, three clones; outputs, decisions and provenance", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-cli-"));
  const files = { "apoe_WT.700.ab1": "apoe_control.ab1", "apoe_A01_hom.700.ab1": "apoe_snp_clone_hom_edited.ab1", "apoe_A02_het.700.ab1": "apoe_snp_clone_het_edited.ab1", "apoe_A03_partial.700.ab1": "apoe_snp_clone_partial_het_edited.ab1" };
  Object.entries(files).forEach(([to, from]) => copyFileSync(sim(from), path.join(dir, to)));
  const out = path.join(dir, "out");
  assert.equal(await cli(["run", dir, "--design", path.join(root, "fixtures", "design", "apoe_r176c_snp.design.json"), "--threads", "2", "--out", out]), 0);
  for (const f of ["results.csv", "results.json", "report.html", "plate.html", "run.json"]) assert.ok(existsSync(path.join(out, f)), `${f} missing`);
  const results = JSON.parse(readFileSync(path.join(out, "results.json"), "utf8")); assert.equal(results.length, 3);
  const by = Object.fromEntries(results.map((r) => [r.name, r]));
  assert.equal(by["apoe_A01_hom.700"].decision.decision, "accept"); assert.equal(by["apoe_A01_hom.700"].well, "A1");
  assert.equal(by["apoe_A02_het.700"].decision.decision, "hold"); assert.notEqual(by["apoe_A03_partial.700"].decision.decision, "accept");
  assert.equal(by["apoe_A01_hom.700"].inputs.edited[0].sha256.length, 64);
  const run_json = JSON.parse(readFileSync(path.join(out, "run.json"), "utf8")); assert.equal(run_json.samples, 3); assert.ok(run_json.decisions.accept >= 1);
  const csv = readFileSync(path.join(out, "results.csv"), "utf8").split("\n")[0]; assert.ok(/decision/.test(csv) && /edited_lo/.test(csv) && /confidence/.test(csv));
  assert.ok(/plate/.test(readFileSync(path.join(out, "plate.html"), "utf8")));
});

test("command line: a sample sheet pairs the files, and a missing file is reported", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-sheet-")); copyFileSync(sim("apoe_control.ab1"), path.join(dir, "ctrlfile.ab1")); copyFileSync(sim("apoe_snp_clone_hom_edited.ab1"), path.join(dir, "c1.ab1"));
  writeFileSync(path.join(dir, "sheet.csv"), "sample,edited,control,type\nclone1,c1,ctrlfile,clone\nclone2,nothere,ctrlfile,clone\n");
  const out = path.join(dir, "out"); const logged = []; const original = console.log; console.log = (...a) => logged.push(a.join(" "));
  try { assert.equal(await cli(["run", dir, "--sheet", path.join(dir, "sheet.csv"), "--design", path.join(root, "fixtures", "design", "apoe_r176c_snp.design.json"), "--threads", "1", "--out", out]), 0); } finally { console.log = original; }
  assert.match(logged.join("\n"), /nothere/);
  assert.equal(JSON.parse(readFileSync(path.join(out, "results.json"), "utf8")).length, 1);
});

test("decisions: the reasons, the thresholds and the cases that must never be accepted", () => {
  const base = { ok: true, kind: "trace-pair", warnings: [], sampleType: "clone", reads: null, agreement: null };
  const iv = (e, i = [0, 0]) => ({ wtPct: [100 - e[1], 100 - e[0]], editedPct: e, intendedEditPct: i, koScorePct: e });
  const homEdit = { ...base, summary: { editedPct: 100, intendedEditPct: 99, wtPct: 0 }, intervals: iv([96, 100], [93, 100]), confidence: { tier: "high", width: 7 }, genotype: { category: "homozygous_edit", summary: "Homozygous for the intended edit", alleles: [] } };
  assert.equal(decide("snp", homEdit).decision, "accept");
  assert.equal(decide("snp", { ...homEdit, confidence: { tier: "low", width: 40 }, intervals: iv([60, 100], [0, 100]) }).decision, "re-sequence");
  assert.equal(decide("snp", { ...homEdit, genotype: { category: "homozygous_partial", summary: "partial", alleles: [] }, summary: { editedPct: 100, intendedEditPct: 1, wtPct: 0 }, intervals: iv([96, 100], [0, 4]) }).decision, "reject");
  assert.equal(decide("snp", { ...homEdit, genotype: { category: "heterozygous_edit", summary: "het", alleles: [] }, intervals: iv([46, 54], [46, 54]) }).decision, "hold");
  assert.equal(decide("snp", { ...homEdit, genotype: { category: "heterozygous_edit", summary: "het", alleles: [] }, intervals: iv([46, 54], [46, 54]) }, { goal: "any" }).decision, "accept");
  assert.equal(decide("snp", { ...homEdit, agreement: { agree: false }, reads: [{ name: "F", editedPct: 90 }, { name: "R", editedPct: 40 }] }).decision, "re-sequence");
  assert.equal(decide("knockout", { ...base, summary: { editedPct: 2, wtPct: 98, intendedEditPct: 0 }, intervals: iv([0, 4]), confidence: { tier: "high", width: 4 }, genotype: { category: "homozygous_wt", summary: "wt", alleles: [] } }).decision, "reject");
  const koHom = { ...base, summary: { editedPct: 100, wtPct: 0, intendedEditPct: 0 }, intervals: iv([95, 100]), confidence: { tier: "high", width: 5 }, genotype: { category: "homozygous_indel", summary: "Homozygous for -1 bp", alleles: [{ kind: "indel", size: -1 }] } };
  assert.equal(decide("knockout", koHom).decision, "accept");
  // a deletion between the cuts that keeps the frame and is shorter than 21 bases is not a knockout (held-out finding), a long one is
  const between = (size) => ({ ...koHom, genotype: { category: "homozygous_deletion", summary: `Homozygous for deletion between cuts (${Math.abs(size)} bp)`, alleles: [{ kind: "deletion_between_cuts", size }] } });
  assert.equal(decide("deletion", between(-18)).decision, "hold");
  assert.equal(decide("deletion", between(-17)).decision, "accept");
  assert.equal(decide("deletion", between(-45)).decision, "accept");
  assert.equal(decide("knockout", { ...koHom, genotype: { ...koHom.genotype, alleles: [{ kind: "indel", size: -3 }] } }).decision, "hold", "an in-frame allele may still make protein");
  assert.equal(decide("knockout", { ok: false, error: "no file" }).decision, "re-sequence");
  assert.deepEqual(tally([{ decision: "accept" }, { decision: "reject" }, { decision: "accept" }]), { accept: 2, hold: 0, "re-sequence": 0, reject: 1 });
});

test("wells are read from names; the plate report lists what to sequence again", () => {
  assert.equal(wellOf("plate1_A01_clone.ab1"), "A1"); assert.equal(wellOf("x-H12.ab1"), "H12"); assert.equal(wellOf("clone14.ab1"), null);
  const results = [{ name: "a", well: "A1", decision: { decision: "accept", reasons: ["ok"], next: "x" } }, { name: "b", well: "B2", decision: { decision: "re-sequence", reasons: ["noisy"], next: "again" } }];
  const html = buildPlateReportHtml({ results }); assert.match(html, /Sequence again \(1\)/); assert.ok(!/<script/i.test(html)); assert.match(resequenceListCsv(results), /"b","B2","noisy"/);
});

test("command line: the figure command writes a self-contained SVG", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-fig-")); const target = path.join(dir, "fig.svg");
  const code = await cli(["figure", "--control", sim("apoe_control.ab1"), "--edited", sim("apoe_snp_clone_hom_edited.ab1"), "--design", path.join(root, "fixtures", "design", "apoe_r176c_snp.design.json"), "--out", target]);
  assert.equal(code, 0); assert.match(readFileSync(target, "utf8"), /^<svg /);
});
