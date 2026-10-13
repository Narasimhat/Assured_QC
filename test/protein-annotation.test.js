// Protein consequence of the designed changes: the reading frame is found from the codons the design names, every designed change gets its amino-acid change,
// and the figures translate the control and each clone at the designed codons. No estimate, call or decision depends on it.
import test from "node:test";
import assert from "node:assert/strict";
import { parseProteinNotes, annotateSpecProtein, notesFromSpec, aaOf, codonsOf } from "../src/protein.js";
import { reverseSpec } from "../src/designSpec.js";
import { parseDesignReport } from "../src/designReport.js";
import { planDesign, analyseSample } from "../src/batch.js";
import { cloneAlignedFigure, clonesAlignedFigure, clonesAlignmentFigure, alignedItem } from "../src/alignedFigure.js";
import { loadSpec, loadTruth, loadTrace, reverseTrace } from "./helpers.js";

const spec = loadSpec("apoe_r176c_snp"); const truth = loadTruth(); const entry = truth.apoe_snp_clone_hom;
const labels = (s) => s.markers.map((m) => (m.protein ? m.protein.label : null));

test("aaOf follows the genetic code", () => {
  assert.equal(aaOf("CGC"), "R"); assert.equal(aaOf("TGC"), "C"); assert.equal(aaOf("TAG"), "*"); assert.equal(aaOf("ATG"), "M"); assert.equal(aaOf("NGC"), "X"); assert.equal(aaOf("AC"), "X");
});

test("notes are read from a design report's text and from a design file's labels", () => {
  const fromReport = parseProteinNotes(["Design: APOE p.R176C", "Silent mutation: p.L177L: CTG -> CTA | PAM TGG->TAG p.D172D: GAC -&gt; GAT | PAM AGG->AGA"]);
  assert.deepEqual(fromReport.map((n) => `${n.from}${n.residue}${n.to} ${n.codonFrom}>${n.codonTo}`), ["L177L CTG>CTA", "D172D GAC>GAT"]);
  const fromSpec = notesFromSpec(spec).map((n) => `${n.from}${n.residue}${n.to} ${n.codonFrom}>${n.codonTo}`);
  assert.deepEqual(fromSpec, ["R176C CGC>TGC", "L173L CTG>CTC", "Q174Q CAG>CAA", "A178A GCA>GCT", "V179V GTG>GTC"]);
  assert.deepEqual(parseProteinNotes(["no codons here", "p.R176C"]), []);
});

test("the frame and strand come from the named codons, and the intended edit gets its residue too", () => {
  const a = annotateSpecProtein(spec);
  assert.equal(a.protein.strand, "+"); assert.deepEqual(labels(a), ["L173L", "Q174Q", "R176C", "A178A", "V179V"]);
  assert.deepEqual(a.markers.map((m) => m.protein.silent), [true, true, false, true, true]);
  assert.equal(a.markers[2].protein.codonFrom, "CGC"); assert.equal(a.markers[2].protein.codonTo, "TGC");
  assert.equal(spec.markers[0].protein, undefined, "the design given is not changed");
  const codons = codonsOf(a, 520, 560); assert.ok(codons.length >= 10); assert.ok(codons.every((c) => c.cols.length === 3 && c.aa.length === 1));
});

test("the same design seen from the other strand gets the same residues", () => {
  const mirrored = reverseSpec(annotateSpecProtein(spec)); const b = annotateSpecProtein(mirrored);
  assert.equal(b.protein.strand, "-"); assert.deepEqual(labels(b).slice().sort(), labels(annotateSpecProtein(spec)).slice().sort());
  assert.deepEqual(b.markers.map((m) => m.protein.codonTo).sort(), ["CAA", "CTC", "GCT", "GTC", "TGC"]);
});

test("without codons in the design, or with codons that fit no frame, the design is returned as it is", () => {
  const bare = { ...spec, design: { ...spec.design, label: "R176C" }, markers: spec.markers.map((m) => ({ ...m, label: "", change: `${m.ref}>${m.alt}` })) };
  assert.equal(annotateSpecProtein(bare), bare);
  const wrong = { ...spec, markers: spec.markers.map((m) => ({ ...m, label: "p.L173L", change: "AAA>AAG" })), design: { ...spec.design, label: "" } };
  assert.equal(annotateSpecProtein(wrong), wrong);
  assert.equal(annotateSpecProtein(null), null);
});

test("a design report carries the codon notes of its silent changes", () => {
  const letters = (s) => [...s].map((b) => `<span>${b}</span>`).join(""); const donor = "ACGT".repeat(15);
  const html = `<h1>Design: APOE p.R176C</h1><h2>2. gRNA Sequences</h2><table><tr><th>Name</th><th>Sequence</th></tr><tr><td>g1</td><td><span>ACGTACGTACGTACGTACGA</span> <span>TGG</span></td></tr></table><h2>4. ssODN Donor Templates</h2><h3>ssODN1</h3><p class="note"><strong>Silent mutation:</strong><br/>p.L177L: CTG -&gt; CTA | PAM TGG-&gt;TAG</p><div><div>Order this strand</div><div><div>Donor</div><div>${letters(donor)}</div></div></div><h2>5. Matched</h2>`;
  const [design] = parseDesignReport(html);
  assert.deepEqual(design.protein.map((n) => `${n.from}${n.residue}${n.to} ${n.codonFrom}>${n.codonTo}`), ["L177L CTG>CTA"]);
});

test("planning a design file annotates it; a manual design takes the notes of a report", () => {
  const planned = planDesign({ designSpec: spec }); assert.deepEqual(labels(planned.spec), ["L173L", "Q174Q", "R176C", "A178A", "V179V"]);
  const control = loadTrace(entry.control); const notes = parseProteinNotes(["p.L173L: CTG -> CTC", "p.Q174Q: CAG -> CAA", "p.A178A: GCA -> GCT", "p.V179V: GTG -> GTC", "R176C (CGC>TGC)"]);
  const manual = planDesign({ control, guides: entry.guides.join(","), donor: entry.donor, gene: "APOE", protein: notes });
  assert.equal(manual.error, undefined, manual.error); assert.deepEqual(labels(manual.spec), ["L173L", "Q174Q", "R176C", "A178A", "V179V"]);
  const plain = planDesign({ control, guides: entry.guides.join(","), donor: entry.donor, gene: "APOE" }); assert.ok(plain.spec.markers.every((m) => !m.protein) && !plain.spec.protein);
});

const analyse = (name, specName, controlTrace, editedTrace) => analyseSample({ name, spec: loadSpec(specName), control: controlTrace, edited: editedTrace, sampleType: "clone" });

test("the view carries the codons, the figures draw the protein row and a box under each clone", () => {
  const r = analyse("apoe_hom", "apoe_r176c_snp", loadTrace(entry.control), loadTrace(entry.edited)); const v = r.view;
  assert.ok(r.ok && v.protein && v.protein.codons.length >= 10); assert.ok(v.markers.every((m) => m.protein && m.protein.cols.length === 3));
  assert.deepEqual(r.markers.map((m) => m.protein && m.protein.label), ["L173L", "Q174Q", "R176C", "A178A", "V179V"]);
  assert.deepEqual(JSON.parse(JSON.stringify(v)), v, "plain data");
  const item = { ...alignedItem({ ...r, workflow: "snp", decision: { decision: "accept" } }) };
  for (const svg of [clonesAlignedFigure([item]), cloneAlignedFigure(v, { name: "hom", decision: "accept" }), clonesAlignmentFigure([item])]) {
    assert.match(svg, /R176C CGC&gt;TGC/); assert.match(svg, /L173L silent CTG&gt;CTC/); assert.match(svg, />protein</); assert.match(svg, /Protein row/);
  }
  assert.match(clonesAlignedFigure([item]), />C176</, "the clone has the intended edit: Cys at residue 176");
  const wt = analyse("apoe_wt", "apoe_r176c_snp", loadTrace(entry.control), loadTrace(entry.control));
  assert.match(clonesAlignedFigure([alignedItem({ ...wt, workflow: "snp", decision: { decision: "reject" } })]), />R176</, "an unedited clone keeps Arg");
});

test("a read from the opposite primer gives the same residues", () => {
  const r = analyse("apoe_rev", "apoe_r176c_snp", reverseTrace(loadTrace(entry.control)), reverseTrace(loadTrace(entry.edited)));
  assert.ok(r.ok, r.error); assert.ok(r.view.protein, "protein row present"); assert.deepEqual(r.markers.map((m) => m.protein && m.protein.label).sort(), ["A178A", "L173L", "Q174Q", "R176C", "V179V"]);
});

test("bases past the stop codon, or labelled noncoding, get no amino-acid change", () => {
  const nr = annotateSpecProtein(loadSpec("nr2f2_2xha_ssodn"));
  assert.deepEqual(nr.markers.map((m) => (m.protein ? m.protein.label : null)), ["P409P", null]);
  const noted = annotateSpecProtein({ ...spec, markers: spec.markers.map((m) => (m.pos === 550 ? { ...m, label: "noncoding" } : m)) });
  assert.equal(noted.markers[4].protein, undefined); assert.equal(noted.markers[3].protein.label, "A178A");
});

test("a design without codon notes draws the figures as before", () => {
  const bare = { ...spec, design: { ...spec.design, label: "R176C" }, markers: spec.markers.map((m) => ({ ...m, label: "", change: `${m.ref}>${m.alt}` })) };
  const r = analyseSample({ name: "bare", spec: bare, control: loadTrace(entry.control), edited: loadTrace(entry.edited), sampleType: "clone" });
  assert.ok(r.ok && r.view); assert.equal(r.view.protein, undefined); assert.ok(r.markers.every((m) => m.protein === undefined));
  const svg = clonesAlignedFigure([alignedItem({ ...r, workflow: "snp", decision: { decision: "accept" } })]);
  assert.doesNotMatch(svg, /Protein row/); assert.doesNotMatch(svg, />protein</);
});
