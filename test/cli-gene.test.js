import test from "node:test";
import assert from "node:assert/strict";
import { mkdtempSync, copyFileSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { main as cli } from "../cli/assured-qc.mjs";
import { validateDesignSpec } from "../src/designSpec.js";
import { planDesign } from "../src/batch.js";
import { spec, gb, geneGenomic, codingSegments, windowStartInGene, designRoot as root } from "./helpers/synthetic-gene.js";

const sim = (n) => path.join(root, "fixtures", "sim", n);
const REF_OFFSET = 1000;
// the design file as the design app writes it for a two-guide knockout: the reference window plus the gene region and its coding segments
const koSpec = () => ({ ...JSON.parse(JSON.stringify(spec)), design: { ...spec.design, editKind: "deletion" }, reference: { ...spec.reference, offset: REF_OFFSET },
  coding: { strand: "+", gene: "NR2F2", transcript: "NM_SYN3.1", offset: REF_OFFSET - windowStartInGene, sequence: geneGenomic, segments: codingSegments } });

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
});

test("design file: the coding block is checked against the reference window and the segments", () => {
  assert.equal(validateDesignSpec(koSpec()).ok, true, JSON.stringify(validateDesignSpec(koSpec()).errors));
  const shifted = koSpec(); shifted.coding.offset += 1; assert.match(validateDesignSpec(shifted).errors.join(" "), /does not contain the reference window/);
  const overlap = koSpec(); overlap.coding.segments = [[10, 40], [30, 60]]; assert.match(validateDesignSpec(overlap).errors.join(" "), /in order and must not overlap/);
  const odd = koSpec(); odd.coding.segments = [[10, 41]]; assert.ok(validateDesignSpec(odd).warnings.some((w) => /not a multiple of three/.test(w)));
  const none = koSpec(); delete none.coding; assert.equal(validateDesignSpec(none).ok, true, "the block is optional");
});

test("design file: planDesign uses the coding block for knockout and deletion designs only", () => {
  const trace = null; // planDesign with a design spec does not read the trace
  const plan = planDesign({ designSpec: koSpec(), control: trace, guides: [], donor: "", gene: "x", nuclease: "SpCas9" });
  assert.ok(!plan.error, plan.error); assert.ok(plan.codingInput, "coding context built from the design file"); assert.equal(plan.codingSource, "design file"); assert.ok(plan.warnings.some((w) => /from the design file/.test(w)));
  const snp = koSpec(); snp.design.editKind = "snp"; assert.equal(planDesign({ designSpec: snp, control: trace, guides: [], donor: "", gene: "x", nuclease: "SpCas9" }).codingInput || null, null);
  const withRecord = planDesign({ designSpec: koSpec(), control: trace, guides: [], donor: "", gene: "x", nuclease: "SpCas9", genbank: gb }); assert.equal(withRecord.codingSource, "gene record", "an uploaded record wins over the design file");
});

test("command line: a design file with a coding block gives the functional assessment without a gene record", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-design-coding-"));
  copyFileSync(sim("nr2f2_control.ab1"), path.join(dir, "NR2F2_WT.700.ab1")); copyFileSync(sim("nr2f2_ko_clone_hom_m1_edited.ab1"), path.join(dir, "NR2F2_A01_hom.700.ab1"));
  const designPath = path.join(dir, "design.json"); writeFileSync(designPath, JSON.stringify(koSpec())); const out = path.join(dir, "out");
  assert.equal(await cli(["run", dir, "--design", designPath, "--workflow", "knockout", "--threads", "1", "--out", out]), 0);
  const results = JSON.parse(readFileSync(path.join(out, "results.json"), "utf8")); assert.equal(results.length, 1);
  assert.ok(results[0].functional, "functional assessment from the design file"); assert.equal(results[0].functional.call.status, "predicted_knockout");
  const run = JSON.parse(readFileSync(path.join(out, "run.json"), "utf8")); assert.equal(run.geneRecord.file, "coding block of the design file"); assert.match(run.geneRecord.sha256, /^[0-9a-f]{64}$/);
  const svg = readFileSync(path.join(out, "protein_products.svg"), "utf8"); assert.match(svg, /^<svg /); assert.match(svg, /Full-length protein/); assert.match(svg, /aa, about [\d.]+ kDa/); assert.ok(!/<script/i.test(svg), "the run figure is written next to the report");
});

test("design file written by the design app (with a coding block) is accepted and gives a coding context", () => {
  const fromApp = JSON.parse(readFileSync(path.join(root, "fixtures", "design", "tagme_ko_deletion_coding.design.json"), "utf8"));
  assert.equal(validateDesignSpec(fromApp).ok, true);
  const plan = planDesign({ designSpec: fromApp, control: null, guides: [], donor: "", gene: "x", nuclease: "SpCas9" });
  assert.equal(plan.codingSource, "design file"); assert.equal(plan.coding.codingExons, 2); assert.equal(plan.coding.proteinLength, 399);
});

test("command line: calibrate scores the functional calls of a run against a protein table", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-calibrate-"));
  copyFileSync(sim("nr2f2_control.ab1"), path.join(dir, "NR2F2_WT.700.ab1")); copyFileSync(sim("nr2f2_ko_clone_hom_m1_edited.ab1"), path.join(dir, "NR2F2_A01_hom.700.ab1")); copyFileSync(sim("nr2f2_ko_clone_het_wt_m1_edited.ab1"), path.join(dir, "NR2F2_A02_het.700.ab1"));
  const designPath = path.join(dir, "design.json"); writeFileSync(designPath, JSON.stringify(koSpec())); const out = path.join(dir, "out");
  assert.equal(await cli(["run", dir, "--design", designPath, "--workflow", "knockout", "--threads", "1", "--out", out]), 0);
  const obs = path.join(dir, "obs.csv"); writeFileSync(obs, "sample,protein,rna\nNR2F2_A01_hom.700,absent,decayed\nNR2F2_A02_het.700,present,\n");
  const report = path.join(dir, "calibration.md"); assert.equal(await cli(["calibrate", path.join(out, "results.json"), obs, "--out", report]), 0);
  const text = readFileSync(report, "utf8"); assert.match(text, /2 clones with a protein result/); assert.match(text, /predicted_knockout/);
  const parsed = JSON.parse(readFileSync(path.join(dir, "calibration.json"), "utf8")); assert.equal(parsed.knockoutCall.called, 1); assert.equal(parsed.knockoutCall.confirmed, 1);
  // the same table decides the clones: protein absent accepts the knockout, protein present rejects the clone that carries wild type
  const rec = path.join(dir, "rec.csv"); assert.equal(await cli(["recommend", path.join(out, "results.json"), "--observations", obs, "--out", rec]), 0);
  const recJson = JSON.parse(readFileSync(path.join(dir, "rec.json"), "utf8")); const by = Object.fromEntries(recJson.map((r) => [r.name, r]));
  assert.equal(by["NR2F2_A01_hom.700"].decision.decision, "accept"); assert.equal(by["NR2F2_A01_hom.700"].recommendation.priority, "confirmed"); assert.equal(by["NR2F2_A02_het.700"].decision.decision, "reject");
});
