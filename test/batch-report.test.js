import test from "node:test";
import assert from "node:assert/strict";
import { analyseSample, analyseJunction, headline, planDesign, bandsFor, WORKFLOWS } from "../src/batch.js";
import { buildQcReportHtml, buildResultsCsv } from "../src/report.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const truth = loadTruth();
const entry = truth.apoe_snp_clone_partial_het;
const spec = loadSpec("apoe_r176c_snp");
const sample = () => analyseSample({ name: "clone <A1>", spec, control: loadTrace(entry.control), edited: loadTrace(entry.edited), sampleType: "clone" });

test("workflows name the method each uses", () => {
  assert.equal(WORKFLOWS.reporter.method, "junction");
  assert.ok(["knockout", "deletion", "snp", "tag_small"].every((key) => WORKFLOWS[key].method === "trace"));
});

test("a design file is validated; a broken one is refused with the first problem", () => {
  assert.equal(planDesign({ designSpec: spec }).source, "design");
  const broken = JSON.parse(JSON.stringify(spec)); broken.markers[0].alt = broken.markers[0].ref;
  assert.match(planDesign({ designSpec: broken }).error, /not usable: Marker at/);
  assert.match(planDesign({}).error, /control trace/);
  const manual = planDesign({ control: loadTrace(entry.control), guides: entry.guides.join(","), donor: entry.donor, gene: "APOE" });
  assert.equal(manual.source, "manual");
  assert.match(planDesign({ control: loadTrace(entry.control), guides: "ACGTACGTACGTACGTACGT" }).error, /not found in the control/);
});

test("a sample result is plain data with charts and a one-line headline", () => {
  const r = sample();
  assert.ok(r.ok);
  assert.equal(JSON.parse(JSON.stringify(r)).genotype.category, "edit_plus_partial");
  assert.match(r.charts.contributions, /^<svg/);
  assert.match(r.charts.discordance, /^<svg/);
  assert.match(headline("snp", r), /partial conversion/);
  const pool = analyseSample({ name: "pool", spec, control: loadTrace(entry.control), edited: loadTrace(truth.apoe_snp_hdr_pool.edited), sampleType: "pool" });
  assert.match(headline("snp", pool), /^Intended edit \d+(\.\d+)?%, blocking-only/);
  assert.equal(pool.genotype, null);
});

test("failures are reported, not thrown", () => {
  const bad = analyseSample({ name: "wrong", spec: loadSpec("nr2f2_sd40v5_block"), control: loadTrace(entry.control), edited: loadTrace(entry.edited) });
  assert.equal(bad.ok, false);
  assert.equal(headline("snp", bad), bad.error);
});

test("report: escapes sample names, states the limits, includes every sample, flags review", () => {
  const html = buildQcReportHtml({ workflow: "snp", design: "APOE R176C", results: [sample(), { name: "x", kind: "trace-pair", ok: false, error: "no <good> trace", warnings: [] }], generatedAt: new Date("2026-10-08T00:00:00Z") });
  assert.ok(html.includes("clone &lt;A1&gt;"));
  assert.ok(!html.includes("clone <A1>"));
  assert.ok(html.includes("no &lt;good&gt; trace"));
  assert.ok(html.includes("Review required"));
  assert.ok(html.includes("Detection limit is about 5%"));
  assert.ok(html.includes("2026-10-08"));
  assert.ok(html.includes("Per-position readout of the designed changes"));
  assert.ok((html.match(/<svg/g) || []).length >= 2);
});

test("csv: one row per sample, quoted fields, failures included", () => {
  const csv = buildResultsCsv("snp", [sample(), { name: "bad, name", ok: false, error: 'bad "trace"', warnings: [] }]);
  const lines = csv.trim().split("\n");
  assert.equal(lines.length, 3);
  assert.ok(lines[0].startsWith("sample,well,decision,decision_reason,status,result"));
  assert.ok(lines[1].includes("edit_plus_partial"));
  assert.ok(lines[2].startsWith('"bad, name",,,,failed'));
  assert.ok(lines[2].includes('""trace""'));
});

test("junction read result and band classification through the batch layer", () => {
  const block = loadSpec("nr2f2_sd40v5_block");
  const r = analyseJunction({ name: "5' read", spec: block, trace: loadTrace("junction5_ok.ab1") });
  assert.ok(r.ok);
  assert.match(headline("reporter", r), /^Matches the expected knock-in allele/);
  const html = buildQcReportHtml({ workflow: "reporter", results: [r] });
  assert.ok(html.includes("5' junction spanned: yes"));
  const bands = bandsFor(block, [841, 631]);
  assert.equal(bands.call, "knock-in and wild type");
  assert.equal(bandsFor({ amplicons: [] }, [1]), null);
});
