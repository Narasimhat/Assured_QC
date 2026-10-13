import test from "node:test";
import assert from "node:assert/strict";
import { specFromControl } from "../src/manualSpec.js";
import { planDesign, analyseSample } from "../src/batch.js";
import { reportDonors, planTasks } from "../src/quick.js";
import { cloneSelectionFigure } from "../src/selectionFigure.js";
import { selectionSection, buildQcReportHtml } from "../src/report.js";
import { withRecommendations, panelOf } from "../src/recommend.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const truth = loadTruth(); const apoe = loadSpec("apoe_r176c_snp");
// two ssODNs for one SNP knock-in: both carry the edit, each carries the blocking change of its own guide
const ref = apoe.reference.sequence; const edit = apoe.markers.find((m) => m.role === "intended"); const blocking = apoe.markers.filter((m) => m.role === "blocking");
const donorWith = (positions, from, to) => { const s = ref.split(""); positions.forEach((p) => { s[p] = apoe.markers.find((m) => m.pos === p).alt; }); return s.slice(from, to).join(""); };
const center = edit.pos; const A = donorWith([edit.pos, blocking[0].pos], center - 60, center + 60); const B = donorWith([edit.pos, blocking[blocking.length - 1].pos], center - 60, center + 60);
const guides = apoe.guides.map((g) => g.spacer);
const control = loadTrace("apoe_control.ab1");

test("several donors are used together: markers merged, the shared change is the intended edit, each donor lists its own markers", () => {
  const built = specFromControl({ control, guides, donor: [{ name: "ssODN A", sequence: A }, { name: "ssODN B", sequence: B }] });
  assert.ok(!built.error, built.error); const spec = built.spec;
  assert.equal(spec.donors.length, 2); assert.deepEqual(spec.donors.map((d) => d.name), ["ssODN A", "ssODN B"]);
  const roles = Object.fromEntries(spec.markers.map((m) => [m.change + m.pos, m.role])); const intended = spec.markers.filter((m) => m.role === "intended");
  assert.equal(intended.length, 1, "only the change both donors carry"); assert.equal(spec.markers.length, 3);
  spec.donors.forEach((d) => { assert.ok(Array.isArray(d.carriesMarkers) && d.carriesMarkers.length === 2); assert.ok(d.carriesMarkers.includes(intended[0].pos)); });
  assert.notDeepEqual(spec.donors[0].carriesMarkers, spec.donors[1].carriesMarkers); assert.equal(spec.design.editKind, "snp"); assert.ok(roles);
});
test("one donor, as a string or a one-item list, builds the same spec as before", () => {
  const single = specFromControl({ control, guides, donor: A }); const listed = specFromControl({ control, guides, donor: [{ name: "x", sequence: A }] });
  assert.deepEqual(listed.spec, single.spec); assert.equal(single.spec.donors[0].name, "donor"); assert.ok(!("carriesMarkers" in single.spec.donors[0]));
});
test("donors that disagree at a position are refused, and an insert falls back to the first donor", () => {
  const k = edit.pos - (center - 60); const third = ref[edit.pos] === "G" || edit.alt === "G" ? "A" : "G"; const clashing = `${B.slice(0, k)}${third}${B.slice(k + 1)}`;
  const clash = specFromControl({ control, guides, donor: [A, clashing] }); assert.match(clash.error || "", /disagree/);
  const withInsert = A.slice(0, 60) + "GCTAGCTAGCTAGC" + A.slice(60); const both = specFromControl({ control, guides, donor: [A, withInsert] });
  assert.ok(both.spec && both.spec.donors.length === 1 && both.warnings.some((w) => /only the first donor/.test(w)));
});
test("a design report's recommended donors are used together; one stays a plain sequence", () => {
  assert.equal(reportDonors([{ name: "a", sequence: "AAAA", recommended: false }, { name: "b", sequence: "CCCC", recommended: true }]), "CCCC");
  assert.deepEqual(reportDonors([{ name: "a", sequence: "AAAA", recommended: true }, { name: "b", sequence: "CCCC", recommended: true }]), [{ name: "a", sequence: "AAAA" }, { name: "b", sequence: "CCCC" }]);
  assert.equal(reportDonors([{ name: "a", sequence: "AAAA", recommended: false }, { name: "b", sequence: "CCCC", recommended: false }]), "AAAA"); assert.equal(reportDonors([]), "");
  const tasks = planTasks({ pairs: [{ sample: "s", control: ["c.ab1"], edited: ["e.ab1"], confidence: "sheet", reason: "" }], getTrace: () => control, design: { kind: "report", guides, donors: [{ name: "A", sequence: A, recommended: true }, { name: "B", sequence: B, recommended: true }] } });
  assert.equal(tasks.problems.length, 0); assert.equal(tasks.tasks[0].spec.donors.length, 2); assert.equal(tasks.tasks[0].workflow, "snp");
});
test("a clone made with either donor is read as the intended edit when both donors are in the design", () => {
  const plan = planDesign({ control, guides, donor: [{ name: "A", sequence: A }, { name: "B", sequence: B }], gene: "sample" }); assert.ok(!plan.error, plan.error);
  const t = truth.apoe_snp_clone_hom; const r = analyseSample({ name: "hom", spec: plan.spec, control: loadTrace(t.control), edited: loadTrace(t.edited), sampleType: "clone" });
  assert.ok(r.ok); assert.ok(r.summary.intendedEditPct > 80, `intended edit ${r.summary.intendedEditPct}%`);
});

test("default panel: six for an edit (replicates), three for a knockout; an explicit number wins", () => {
  const mk = (name, workflow, lo) => ({ ok: true, kind: "trace-pair", name, workflow, sampleType: "clone", confidence: { tier: "high" }, warnings: [], summary: { wtPct: 0, editedPct: 100, intendedEditPct: 100, koScorePct: 100 }, intervals: { editedPct: [96, 100], intendedEditPct: [lo, 100] },
    genotype: { category: "homozygous_edit", alleles: [{ label: "Intended edit", kind: "edit", size: 0 }] }, contributions: [{ label: "Intended edit", kind: "edit", size: 0, fraction: 1 }], decision: { decision: "accept", reasons: [], next: "" } });
  const edits = Array.from({ length: 8 }, (_, i) => mk(`e${i}`, "snp", 99 - i)); assert.equal(panelOf(withRecommendations(edits)).length, 6); assert.equal(panelOf(withRecommendations(edits, { panel: 4 })).length, 4);
});

test("the clone selection figure and the report section: chromatograms of the clones to carry forward, a table of all clones, next steps", () => {
  const plan = planDesign({ control, guides, donor: A, gene: "sample" });
  const make = (key, name, decision) => { const t = truth[key]; const r = analyseSample({ name, spec: plan.spec, control: loadTrace(t.control), edited: loadTrace(t.edited), sampleType: "clone" }); return { ...r, workflow: "snp", decision: { decision, reasons: [], next: "" } }; };
  const rs = withRecommendations([make("apoe_snp_clone_hom", "clone_hom", "accept"), make("apoe_snp_clone_het", "clone_het", "hold")]);
  const svg = cloneSelectionFigure(rs); assert.match(svg, /^<svg /); assert.ok(!/<script|onload|href=/i.test(svg)); assert.match(svg, /of 2 clones? to carry forward first/); assert.match(svg, />clone_hom</); assert.match(svg, /All clones/); assert.match(svg, /Carry forward first \(/);
  assert.ok(svg.indexOf(">clone_hom<") < svg.lastIndexOf(">clone_het<")); assert.ok((svg.match(/<polyline /g) || []).length >= 8);
  const section = selectionSection(rs, "snp"); assert.match(section, /<h2>Clone selection<\/h2>/); assert.match(section, /Next steps for/); assert.match(section, /Karyotype/);
  const html = buildQcReportHtml({ workflow: "snp", results: rs }); assert.ok(html.indexOf("Clone selection") < html.indexOf("All clones aligned"));
  assert.equal(cloneSelectionFigure([]), ""); assert.equal(selectionSection([{ ok: false, name: "x" }], "snp"), "");
});

test("SNP and small-tag runs name the intended-edit share the HDR score (ICE's KI score): table, figure, headline and ranking basis", () => {
  const plan = planDesign({ control, guides, donor: A, gene: "sample" });
  const make = (key, name, decision) => { const t = truth[key]; const r = analyseSample({ name, spec: plan.spec, control: loadTrace(t.control), edited: loadTrace(t.edited), sampleType: "clone" }); return { ...r, workflow: "snp", decision: { decision, reasons: [], next: "" } }; };
  const rs = withRecommendations([make("apoe_snp_clone_hom", "clone_hom", "accept"), make("apoe_snp_clone_het", "clone_het", "hold")]);
  assert.match(rs[0].recommendation.basis, /^HDR score/);
  const svg = cloneSelectionFigure(rs); assert.match(svg, /HDR score \(95% interval\)/); assert.match(svg, /KI score/);
  const shown = (svg.match(/\d+\.\d \((?:\d+(?:\.\d+)?)\u2013(?:\d+(?:\.\d+)?)\)/g) || []); assert.ok(shown.length >= 2, `a score with its interval for each clone: ${shown}`);
  const html = buildQcReportHtml({ workflow: "snp", results: rs }); assert.match(html, /<th>HDR score \(intended edit\)<\/th>/); assert.match(html, /HDR score \d+(\.\d+)?-\d+(\.\d+)?%/);
  assert.ok(!/KO score/.test(svg), "an edit run does not talk about a KO score");
});
