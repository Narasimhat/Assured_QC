import test from "node:test";
import assert from "node:assert/strict";
import { analysePair } from "../src/analyze.js";
import { specFromControl, locateGuide } from "../src/manualSpec.js";
import { validateDesignSpec } from "../src/designSpec.js";
import { callGenotype } from "../src/genotype.js";
import { loadSpec, loadTruth, loadTrace, reverseTrace } from "./helpers.js";

const truth = loadTruth();
const specFor = (name) => loadSpec(name.startsWith("nr2f2") ? "nr2f2_2xha_ssodn" : "apoe_r176c_snp");
const analyse = (name, extra = {}, tweak = (x) => x) => {
  const entry = truth[name];
  return analysePair({ control: tweak(loadTrace(entry.control)), edited: tweak(loadTrace(entry.edited)), spec: specFor(name), options: { sampleType: "clone", ...extra } });
};
const near = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg}: ${a} vs ${b}`);

test("clone genotype calls on the synthetic clones", () => {
  const expected = {
    nr2f2_ko_clone_het_wt_m1: "heterozygous_indel", nr2f2_ko_clone_compound: "compound_heterozygous", nr2f2_ko_clone_hom_m1: "homozygous_indel",
    nr2f2_ha_hdr_clone_het: "heterozygous_edit", nr2f2_ha_hdr_clone_hom: "homozygous_edit",
    apoe_snp_clone_het: "heterozygous_edit", apoe_snp_clone_hom: "homozygous_edit", apoe_snp_clone_partial_het: "edit_plus_partial",
  };
  for (const [name, category] of Object.entries(expected)) {
    const result = analyse(name);
    assert.ok(result.ok, result.error);
    assert.equal(result.genotype.category, category, `${name}: ${result.genotype.summary}`);
  }
});

test("a pool is not given a genotype", () => {
  const result = analyse("nr2f2_ko_pool_g1", { sampleType: "pool" });
  assert.equal(result.genotype, null);
});

test("marker readout: the SNP and the blocking changes are read independently of the fit", () => {
  const partial = analyse("apoe_snp_clone_partial_het").markerReadout;
  const snp = partial.find((row) => row.marker.role === "intended");
  const silent = partial.filter((row) => row.marker.role === "blocking");
  near(snp.altNet, 0.5, 0.06, "SNP allele fraction");
  silent.forEach((row) => assert.ok(row.altNet > 0.88, `blocking marker ${row.marker.pos}: ${row.altNet}`));
  const het = analyse("apoe_snp_clone_het").markerReadout;
  het.forEach((row) => near(row.altNet, 0.5, 0.08, `het marker ${row.marker.pos}`));
  const hom = analyse("apoe_snp_clone_hom").markerReadout;
  hom.forEach((row) => assert.ok(row.altNet > 0.85, `hom marker ${row.marker.pos}: ${row.altNet}`));
  const wildtype = analyse("nr2f2_ko_clone_het_wt_m1").markerReadout.filter((row) => row.covered);
  wildtype.forEach((row) => assert.ok(row.altNet < 0.03, "no donor base in a knockout clone"));
});

test("markers past an insertion are not read position by position", () => {
  const rows = analyse("nr2f2_ha_hdr_clone_het").markerReadout;
  assert.ok(rows.some((row) => !row.covered && /downstream of the insertion/.test(row.reason)));
  assert.ok(rows.some((row) => row.covered && row.altNet > 0.4 && row.altNet < 0.6));
});

test("a read from the reverse primer gives the same answer as one from the forward primer", () => {
  const forward = analyse("nr2f2_ko_pool_g1", { sampleType: "pool" });
  const reverse = analyse("nr2f2_rev_ko_pool_g1", { sampleType: "pool" });
  assert.ok(reverse.ok, reverse.error);
  assert.equal(forward.prepared.orientation, "+");
  assert.equal(reverse.prepared.orientation, "-");
  near(reverse.summary.editedPct, 75, 4, "edited (reverse read)");
  near(reverse.summary.wtPct, 25, 4, "wild type (reverse read)");
  near(reverse.summary.koScorePct, 53, 4, "KO score (reverse read)");
  near(reverse.summary.editedPct, forward.summary.editedPct, 4, "forward and reverse reads agree");
  assert.equal(reverse.prepared.spec.guides[0].strand, "-", "the design was turned around, not the trace");
  // a trace turned around sample by sample is the same data read the other way
  const turned = analyse("apoe_snp_clone_partial_het", {}, reverseTrace);
  assert.ok(turned.ok, turned.error);
  assert.equal(turned.prepared.orientation, "-");
  assert.equal(turned.genotype.category, "edit_plus_partial");
});

test("a reference that differs from the control near the cut is patched, with a warning", () => {
  const spec = specFor("nr2f2_ko_pool_g1");
  const entry = truth.nr2f2_ko_pool_g1;
  const start = spec.reference.sequence.indexOf(entry.reference);
  const position = start + entry.cuts.g1 - 12; // 12 bases upstream of the cut
  const original = spec.reference.sequence[position];
  const changed = original === "A" ? "C" : "A";
  const altered = { ...spec, reference: { ...spec.reference, sequence: spec.reference.sequence.slice(0, position) + changed + spec.reference.sequence.slice(position + 1) } };
  const result = analysePair({ control: loadTrace(entry.control), edited: loadTrace(entry.edited), spec: altered, options: { sampleType: "pool" } });
  assert.ok(result.ok, result.error);
  assert.ok(result.prepared.patches.length >= 1);
  assert.ok(result.warnings.some((w) => /control differs from the design reference/.test(w)));
  near(result.summary.editedPct, 75, 5, "editing is still found");
});

test("clear failures: wrong design, guide outside the read, edited trace from a different amplicon", () => {
  const entry = truth.apoe_snp_clone_het;
  const wrong = analysePair({ control: loadTrace(entry.control), edited: loadTrace(entry.edited), spec: loadSpec("nr2f2_sd40v5_block") });
  assert.equal(wrong.ok, false);
  assert.match(wrong.error, /does not match the design reference/);
  const mismatched = analysePair({ control: loadTrace("apoe_control.ab1"), edited: loadTrace("nr2f2_ko_pool_g1_edited.ab1"), spec: specFor("apoe_snp_clone_het") });
  assert.equal(mismatched.ok, false);
  assert.match(mismatched.error, /edited trace does not align to the control/);
  const noGuides = analysePair({ control: loadTrace("apoe_control.ab1"), edited: loadTrace(entry.edited), spec: { ...specFor("apoe_snp_clone_het"), guides: [] } });
  assert.match(noGuides.error, /no guides/);
});

test("manual mode: guide and donor entered by hand reproduce the design-based answer", () => {
  const entry = truth.apoe_snp_clone_partial_het;
  const donor = entry.donor;
  const built = specFromControl({ control: loadTrace(entry.control), guides: entry.guides.join(","), donor, gene: "APOE" });
  assert.equal(built.error, undefined, built.error);
  assert.deepEqual(validateDesignSpec(built.spec).errors, []);
  assert.equal(built.spec.markers.length, 5);
  assert.equal(built.spec.design.editKind, "snp");
  const result = analysePair({ control: loadTrace(entry.control), edited: loadTrace(entry.edited), spec: built.spec, options: { sampleType: "clone" } });
  assert.ok(result.ok, result.error);
  // markers are unclassified in manual mode: the complete edit is still separated from the partial one
  near(result.summary.intendedEditPct, 50, 8, "complete edit");
  near(result.summary.partialConversionPct, 50, 8, "partial conversion");
});

test("manual mode input checks", () => {
  const control = loadTrace("nr2f2_control.ab1");
  assert.match(specFromControl({ control, guides: "" }).error, /at least one guide/);
  assert.match(specFromControl({ control, guides: "ACGTACGTACGTACGTACGT" }).error, /not found in the control/);
  assert.match(specFromControl({ control, guides: "ACGT" }).error, /17-24 nt/);
  assert.match(specFromControl({ control, guides: "CAGUUUUAACUGGCCGUAUA", donor: "A".repeat(301) }).error, /longer than 300/);
  const built = specFromControl({ control, guides: "cagttttaactggccgtata" });
  assert.equal(built.spec.guides[0].strand, "+");
  assert.equal(built.spec.guides[0].pam, "TGG");
  assert.equal(locateGuide(control.calls, "AATAAATAAATAAAATAAGA").cut - locateGuide(control.calls, "CAGTTTTAACTGGCCGTATA").cut, 29);
});

test("genotype rules on constructed decompositions", () => {
  const contribution = (kind, fraction, size = 0, extra = {}) => ({ kind, fraction, size, label: kind, ...extra });
  const call = (contributions) => callGenotype({ contributions }).category;
  assert.equal(call([contribution("edit", 0.92)]), "homozygous_edit");
  assert.equal(call([contribution("edit", 0.52), contribution("wt", 0.46)]), "heterozygous_edit");
  assert.equal(call([contribution("edit", 0.5), contribution("edit_partial", 0.45, 0, { label: "Partial conversion: x" })]), "edit_plus_partial");
  assert.equal(call([contribution("indel", 0.48, -1), contribution("indel", 0.47, 4)]), "compound_heterozygous");
  assert.equal(call([contribution("wt", 0.4), contribution("indel", 0.3, -1), contribution("indel", 0.3, 2)]), "mixed");
  assert.equal(call([contribution("wt", 0.97)]), "homozygous_wt");
  assert.equal(call([contribution("wt", 0.004)]), "unclear");
});

test("three guides: the library and the fit handle them, and the deletion between the outer cuts is found", () => {
  const entry = truth.nr2f2_del_two_guides;
  const built = specFromControl({ control: loadTrace(entry.control), guides: [...entry.guides, "ATTTATTGAATTGCCATATA"].join(","), gene: "NR2F2" });
  assert.equal(built.error, undefined, built.error);
  assert.equal(built.spec.guides.length, 3);
  const result = analysePair({ control: loadTrace(entry.control), edited: loadTrace(entry.edited), spec: built.spec, options: { sampleType: "pool" } });
  assert.ok(result.ok, result.error);
  near(result.summary.editedPct, 70, 8, "edited");
  const between = result.decomposition.contributions.find((c) => c.kind === "deletion_between_cuts");
  assert.ok(between && Math.abs(between.size + 29) <= 3, `deletion between cuts: ${between?.size}`);
  assert.match(specFromControl({ control: loadTrace(entry.control), guides: ["CAGTTTTAACTGGCCGTATA", "AATAAATAAATAAAATAAGA", "ATTTATTGAATTGCCATATA", "ACACACCTCATGTGACCCAA"] }).error, /At most three guides/);
});

test("systematic trace artefacts shared by control and edited do not look like editing", () => {
  const specNr = loadSpec("nr2f2_2xha_ssodn");
  const run = (name) => analysePair({ control: loadTrace(truth[name].control), edited: loadTrace(truth[name].edited), spec: specNr, options: { sampleType: "clone" } });
  const nothing = run("nr2f2_artifact_null");
  assert.ok(nothing.ok, nothing.error);
  assert.ok(nothing.summary.editedPct < 2, `spurious editing ${nothing.summary.editedPct}%`);
  assert.equal(nothing.genotype.category, "homozygous_wt");
  const het = run("nr2f2_artifact_ko_het");
  near(het.summary.wtPct, 50, 5, "wild type");
  assert.equal(het.genotype.category, "heterozygous_indel");
  // without the control's own composition for unchanged context, the same artefacts are misread as editing
  const naive = analysePair({ control: loadTrace(truth.nr2f2_artifact_null.control), edited: loadTrace(truth.nr2f2_artifact_null.edited), spec: specNr, options: { controlContext: false } });
  assert.ok(naive.summary.editedPct > 5, "the artefact fixture is strong enough to matter");
});

test("noisier traces: a one- or two-position edit is not pulled down by the noise handling", () => {
  const specApoe = loadSpec("apoe_r176c_snp");
  const run = (name) => analysePair({ control: loadTrace(truth[name].control), edited: loadTrace(truth[name].edited), spec: specApoe, options: { sampleType: "clone" } });
  const het = run("apoe_noisy_snp_het");
  assert.ok(het.ok, het.error);
  near(het.summary.intendedEditPct, 50, 6, "het clone intended edit");
  near(het.summary.wtPct, 50, 6, "het clone wild type");
  assert.equal(het.genotype.category, "heterozygous_edit");
  const hom = run("apoe_noisy_snp_hom");
  assert.ok(hom.summary.intendedEditPct >= 90, `hom intended ${hom.summary.intendedEditPct}`);
  assert.equal(hom.genotype.category, "homozygous_edit");
});
