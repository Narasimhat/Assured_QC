import test from "node:test";
import assert from "node:assert/strict";
import { runSummaryFigure, summaryOrder } from "../src/summaryFigure.js";
import { withRecommendations } from "../src/recommend.js";

const row = (label, fraction, kind = "deletion_between_cuts", size = -79) => ({ label, kind, size, fraction });
const clone = (name, contributions, over = {}) => ({ ok: true, kind: "trace-pair", name, workflow: "deletion", sampleType: "clone", confidence: { tier: "high" }, warnings: [],
  summary: { wtPct: 0, editedPct: 100, koScorePct: 100 }, intervals: { editedPct: [96, 100], wtPct: [0, 3] }, contributions,
  genotype: { category: "homozygous_indel", alleles: contributions.filter((c) => c.kind !== "wt").map((c) => ({ label: c.label, kind: c.kind, size: c.size })) },
  decision: { decision: "hold", reasons: ["caveats"], next: "check protein" }, ...over });
const results = () => [
  clone("b_clone", [row("Deletion between cuts (79 bp)", 0.5), row("Deletion between cuts (77 bp)", 0.5, "deletion_between_cuts", -77)], { warnings: ["The edited trace has low quality scores over much of its length; treat the result with caution."] }),
  clone("a_clone", [row("Deletion between cuts (79 bp)", 0.97)], { decision: { decision: "accept", reasons: ["ok"], next: "" } }),
  clone("c_wt", [{ label: "Wild type", kind: "wt", size: 0, fraction: 0.95 }], { decision: { decision: "reject", reasons: ["wt"], next: "discard" }, summary: { wtPct: 95, editedPct: 5, koScorePct: 0 }, intervals: { editedPct: [1, 9], wtPct: [91, 99] } }),
  { ok: false, name: "d_bad", error: "no usable read", decision: { decision: "re-sequence", reasons: ["no read"], next: "again" } },
];

test("the summary graph has one row per sample, best clone first, with decision, composition and plan", () => {
  const ranked = withRecommendations(results()); const svg = runSummaryFigure(ranked, { title: "Run summary: plate 1" });
  assert.match(svg, /^<svg /); assert.ok(svg.trimEnd().endsWith("</svg>")); assert.match(svg, /Run summary: plate 1/); assert.match(svg, /Accept 1   \u00b7   Hold 1   \u00b7   Sequence again 1   \u00b7   Reject 1/);
  const at = (s) => svg.indexOf(`>${s}</text>`); assert.ok(at("a_clone") > 0 && at("a_clone") < at("b_clone") && at("b_clone") < at("c_wt") && at("c_wt") < at("d_bad"), "ranked clones first, then the rest in input order");
  assert.match(svg, />\u221279</); assert.match(svg, />\u221277</); assert.match(svg, />WT</); assert.match(svg, /#1 carry forward/); assert.match(svg, /#2 carry forward/);
  assert.match(svg, /Deletion between the cuts/); assert.match(svg, /Wild type/); assert.match(svg, /no usable read/); assert.match(svg, /KO score/); assert.ok(!/Knockout support/.test(svg), "no gene record: the column is the KO score");
  assert.equal((svg.match(/<svg /g) || []).length, 1); assert.ok(!/<script|javascript:/i.test(svg));
  assert.match(svg, /<text[^>]*>!<\/text>/, "the low-quality clone is marked");
});

test("the order follows the rank when there is one, and the input order when there is not; the graph still draws without a recommendation", () => {
  const ranked = withRecommendations(results()); const order = summaryOrder(ranked).map((i) => ranked[i].name); assert.deepEqual(order, ["a_clone", "b_clone", "c_wt", "d_bad"]);
  assert.deepEqual(summaryOrder(results()), [0, 1, 2, 3]);
  const plain = runSummaryFigure(results()); assert.match(plain, /b_clone/); assert.doesNotMatch(plain, /carry forward/);
  assert.equal(runSummaryFigure([]), ""); assert.equal(runSummaryFigure(null), "");
});

test("rows shrink for large runs, edit workflows show the intended edit, text is escaped", () => {
  const many = Array.from({ length: 90 }, (_, i) => clone(`clone${i + 1}`, [row("Deletion between cuts (79 bp)", 0.97)]));
  const big = runSummaryFigure(many); const small = runSummaryFigure(many.slice(0, 10)); const h = (s) => Number(/height="(\d+)"/.exec(s)[1]); assert.ok(h(big) < 90 * 28 / 1.5, "rows are 15 px for a run this large"); assert.ok(h(small) > 10 * 27);
  const snp = { ...clone("snp1", [{ label: "Intended edit", kind: "edit", size: 0, fraction: 0.9 }, { label: "Wild type", kind: "wt", size: 0, fraction: 0.1 }]), workflow: "snp", summary: { intendedEditPct: 90, editedPct: 90 }, intervals: { intendedEditPct: [84, 94], editedPct: [84, 94] } };
  const s = runSummaryFigure([snp]); assert.match(s, /HDR score, 95% interval/); assert.match(s, /Lower bound/); assert.match(s, /ICE's KI score/); assert.match(s, /Intended edit/); assert.match(s, />edit</); assert.match(s, /90%/);
  const evil = runSummaryFigure([clone("<img src=x onerror=1>&\"", [row("<b>x</b>", 0.97)])]); assert.ok(!/<img|<b>/.test(evil)); assert.match(evil, /&lt;img src=x onerror=1&gt;&amp;&quot;/);
});
