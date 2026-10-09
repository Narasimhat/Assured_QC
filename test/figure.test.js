import test from "node:test";
import assert from "node:assert/strict";
import { chromatogramFigure } from "../src/figure.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";
const truth = loadTruth();
test("figure: both chromatograms, the cuts and the alleles in one self-contained SVG", () => {
  const entry = truth.nr2f2_ko_clone_het_wt_m1;
  const { svg, result } = chromatogramFigure({ name: "clone 7", controlTrace: loadTrace(entry.control), editedTrace: loadTrace(entry.edited), spec: loadSpec("nr2f2_2xha_ssodn"), sampleType: "clone", workflow: "knockout" });
  assert.match(svg, /^<svg /); assert.ok(!/<script|onload|href=/i.test(svg)); assert.ok(svg.length < 600000, `svg is ${svg.length} bytes`);
  assert.ok((svg.match(/<path /g) || []).length === 8, "four channels in each of two panels");
  assert.match(svg, /clone 7: edited (4|5)\d(\.\d)?%/); assert.match(svg, /95% interval/); assert.match(svg, /Wild type/);
  assert.ok(result.decision && ["accept", "hold", "reject", "re-sequence"].includes(result.decision.decision));
});
test("figure: an unusable pair gives a message, not a broken figure", () => {
  const entry = truth.nr2f2_ko_clone_het_wt_m1; const other = truth.apoe_snp_clone_hom;
  const { svg, result } = chromatogramFigure({ name: "x", controlTrace: loadTrace(entry.control), editedTrace: loadTrace(other.edited), spec: loadSpec("nr2f2_2xha_ssodn") });
  assert.equal(result, null); assert.match(svg, /<text/);
});
