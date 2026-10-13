import test from "node:test";
import assert from "node:assert/strict";
import { analyseSample } from "../src/batch.js";
import { alignWindow } from "../src/alignedView.js";
import { cloneAlignedFigure, clonesAlignedFigure, clonesAlignmentFigure, alignedGroups } from "../src/alignedFigure.js";
import { buildQcReportHtml, alignedSection } from "../src/report.js";
import { withRecommendations } from "../src/recommend.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const truth = loadTruth();
const analyse = (key, specName, options = {}) => { const t = truth[key]; return analyseSample({ name: key, spec: loadSpec(specName), control: loadTrace(t.control), edited: loadTrace(t.edited), sampleType: "clone", options }); };
const decided = (r, decision = "accept") => ({ ...r, workflow: "snp", decision: { decision, reasons: [], next: "" } });
// a fixed pseudo-random sequence, so that every stretch is unique
const seq = (n, seed = 7) => { let s = seed; let out = ""; for (let i = 0; i < n; i += 1) { s = (s * 1103515245 + 12345) % 2147483648; out += "ACGT"[(s >> 16) & 3]; } return out; };

test("alignWindow: identical reads map one to one; a mismatch is not a gap", () => {
  const c = seq(80); assert.deepEqual(alignWindow(c, c, 0, 79), Array.from({ length: 80 }, (_, i) => i));
  const e = `${c.slice(0, 40)}${c[40] === "A" ? "C" : "A"}${c.slice(41)}`; const m = alignWindow(c, e, 0, 79);
  assert.ok(m.every((k, i) => k === i), "a single changed base keeps the mapping");
});
test("alignWindow: one deletion is one gap, one insertion is one jump", () => {
  const c = seq(90); const del = `${c.slice(0, 30)}${c.slice(35)}${seq(20, 3)}`; const m = alignWindow(c, del, 0, del.length - 1);
  assert.deepEqual(m.slice(0, 30), Array.from({ length: 30 }, (_, i) => i)); assert.ok(m.slice(30, 35).every((k) => k === -1), "the five missing bases are the gap");
  assert.deepEqual(m.slice(35), Array.from({ length: 55 }, (_, i) => 30 + i));
  const ins = `${c.slice(0, 50)}GGG${c.slice(50)}${seq(10, 5)}`; const k = alignWindow(c, ins, 0, ins.length - 1);
  assert.ok(k.every((x) => x >= 0)); assert.equal(k.filter((x, i) => i > 0 && x - k[i - 1] > 1).length, 1, "one jump"); assert.ok(k.some((x, i) => i > 0 && x - k[i - 1] === 4), "three extra bases in the clone");
});

test("a SNP knock-in clone: guides with PAM and cut, the designed changes, both reads over the same stretch", () => {
  const r = analyse("apoe_snp_clone_hom", "apoe_r176c_snp"); const v = r.view; const spec = loadSpec("apoe_r176c_snp");
  assert.ok(r.ok && v); assert.equal(v.control.calls.length, v.length); assert.equal(v.edited.calls.length, v.length);
  assert.equal(v.guides.length, spec.guides.length); v.guides.forEach((g) => { assert.equal(g.spacer.length, 2); assert.equal(g.pam.length, 2); assert.equal(g.pam[1] - g.pam[0], 2); });
  assert.equal(v.markers.length, spec.markers.length); assert.equal(v.mode, "offset");
  assert.ok(v.markers.some((m) => v.edited.calls[m.at] !== v.control.calls[m.at]), "the homozygous edit shows as a base that differs from the control");
  for (const b of "ACGT") { assert.equal(v.control.curves[b].length, v.length * v.points); assert.equal(v.edited.curves[b].length, v.length * v.points); assert.ok(v.edited.curves[b].every((x) => x >= 0 && x <= 100)); }
  assert.deepEqual(JSON.parse(JSON.stringify(v)), v, "plain data, so it can cross a worker boundary");
});
test("a clone with a clean deletion on both alleles is drawn with a gap; a mixed read is not forced into one", () => {
  const hom = analyse("nr2f2_ko_clone_hom_m1", "nr2f2_2xha_ssodn"); assert.equal(hom.view.mode, "gapped"); assert.equal((hom.view.edited.calls.match(/-/g) || []).length, 1);
  const het = analyse("nr2f2_ko_clone_het_wt_m1", "nr2f2_2xha_ssodn"); assert.equal(het.view.mode, "offset"); assert.ok(!het.view.edited.calls.includes("-"));
});
test("a read from the reverse primer is drawn in its own direction", () => {
  const r = analyse("nr2f2_rev_ko_pool_g1", "nr2f2_2xha_ssodn"); assert.ok(r.ok); assert.equal(r.view.orientation, "-"); assert.ok(r.view.guides.length >= 1);
});
test("the view can be switched off and never stops an analysis", () => {
  const r = analyse("apoe_snp_clone_hom", "apoe_r176c_snp", { traceView: false }); assert.ok(r.ok); assert.equal(r.view, undefined);
});

test("figures: one clone against the control, all clones stacked, all clones as an alignment", () => {
  const hom = decided(analyse("apoe_snp_clone_hom", "apoe_r176c_snp")); const svg = cloneAlignedFigure(hom.view, { name: "clone 7", decision: "accept", headline: "Homozygous for the intended edit" });
  assert.match(svg, /^<svg /); assert.ok(svg.trimEnd().endsWith("</svg>")); assert.ok(!/<script|onload|href=/i.test(svg)); assert.match(svg, />clone 7</); assert.match(svg, /Wild type/); assert.match(svg, /Accept/);
  assert.match(svg, /gRNA/); assert.ok((svg.match(/<polyline /g) || []).length >= 8, "four channels, control and clone"); assert.ok(svg.length < 150000, `svg is ${svg.length} bytes`);
  const het = { ...decided(analyse("apoe_snp_clone_het", "apoe_r176c_snp"), "hold"), name: "clone_het" }; const ranked = withRecommendations([het, { ...hom, name: "clone_hom" }]); const groups = alignedGroups(ranked);
  assert.equal(groups.length, 1); assert.deepEqual(groups[0].map((x) => x.name), ["clone_hom", "clone_het"], "best clone first");
  const stack = clonesAlignedFigure(groups[0]); const grid = clonesAlignmentFigure(groups[0]);
  for (const s of [stack, grid]) { assert.match(s, /^<svg /); assert.ok(s.indexOf(">clone_hom<") > 0 && s.indexOf(">clone_hom<") < s.indexOf(">clone_het<")); assert.match(s, /Wild type/); }
  assert.ok((grid.match(/<circle /g) || []).length > 100, "a dot for every base that agrees with the control");
});
test("samples of another design form their own group", () => {
  const a = decided(analyse("apoe_snp_clone_hom", "apoe_r176c_snp")); const b = decided(analyse("nr2f2_ko_clone_hom_m1", "nr2f2_2xha_ssodn"));
  assert.equal(alignedGroups(withRecommendations([a, b])).length, 2);
});

test("report: every clone has its aligned chromatogram and the run has an all-clones section", () => {
  const rs = withRecommendations([decided(analyse("apoe_snp_clone_hom", "apoe_r176c_snp")), decided(analyse("apoe_snp_clone_het", "apoe_r176c_snp"), "hold")]);
  const html = buildQcReportHtml({ workflow: "snp", results: rs });
  assert.match(html, /<h2>All clones aligned<\/h2>/); assert.equal((html.match(/Chromatogram aligned with the wild-type control/g) || []).length, 2); assert.match(html, /Base calls/);
  const off = buildQcReportHtml({ workflow: "snp", results: rs.map((r) => ({ ...r, view: undefined })) }); assert.ok(!/All clones aligned|Chromatogram aligned/.test(off));
});
test("report: above 48 clones the stacked chromatograms are left out and wild-type rejects get no per-sample figure", () => {
  const one = decided(analyse("apoe_snp_clone_hom", "apoe_r176c_snp"), "reject"); const many = Array.from({ length: 50 }, (_, i) => ({ ...one, name: `c${i}` }));
  const section = alignedSection(many); assert.match(section, /More than 48 clones/); assert.ok(!/Chromatograms<\/h3>/.test(section));
  const html = buildQcReportHtml({ workflow: "snp", results: many }); assert.ok(!/Chromatogram aligned with the wild-type control/.test(html));
  const kept = buildQcReportHtml({ workflow: "snp", results: [...many.slice(0, 49), { ...one, name: "keep", decision: { decision: "hold", reasons: [], next: "" } }] });
  assert.equal((kept.match(/Chromatogram aligned with the wild-type control/g) || []).length, 1);
});
