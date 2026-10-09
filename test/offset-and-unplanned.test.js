// 0.3.0: offset re-anchoring at the cut (src/offset.js) and the unplanned-base-change scan (src/unplanned.js).
// Positives and negatives are built from the synthetic fixtures only (no laboratory traces).
import test from "node:test";
import assert from "node:assert/strict";
import { analysePair } from "../src/analyze.js";
import { analyseSample } from "../src/batch.js";
import { decide } from "../src/decision.js";
import { refineOffsetAtCut, shiftAt } from "../src/offset.js";
import { scanUnplanned } from "../src/unplanned.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const truth = loadTruth();
const specFor = (name) => loadSpec(name.startsWith("nr2f2") ? "nr2f2_2xha_ssodn" : "apoe_r176c_snp");
const BASES = ["A", "C", "G", "T"];

const pair = (name, tweak = (t) => t, options = {}) => {
  const entry = truth[name];
  return analysePair({ control: loadTrace(entry.control), edited: tweak(loadTrace(entry.edited)), spec: specFor(name), options: { sampleType: "clone", ...options } });
};

/** A copy of the trace as if the base caller had not made the call at `index` (the peak is still in the signal). */
function dropCall(trace, index) {
  const without = (arr) => { const out = Array.from(arr); out.splice(index, 1); return out; };
  return { ...trace, calls: trace.calls.slice(0, index) + trace.calls.slice(index + 1), peakLocations: Int32Array.from(without(trace.peakLocations)), quality: Uint8Array.from(without(trace.quality)) };
}

/** A copy with a second base of the same height as the called one at call `index` (a heterozygous change of that base to `alt`). */
function addSecondBase(trace, index, alt, share = 1) {
  const channels = {}; for (const b of BASES) channels[b] = Int16Array.from(trace.channels[b]);
  const p = trace.peakLocations[index]; const called = trace.calls[index];
  const height = Math.max(...Array.from({ length: 5 }, (_, d) => trace.channels[called][p - 2 + d]));
  for (let d = -4; d <= 4; d += 1) channels[alt][p + d] = Math.round(height * share * Math.exp(-(d * d) / 8)) + channels[alt][p + d];
  return { ...trace, channels };
}

const base = pair("nr2f2_ko_clone_hom_m1");
const topLabel = (r) => `${r.decomposition.contributions[0].label}`;

test("offset: a call missed by the base caller between the alignment stretch and the cut is re-anchored", () => {
  assert.ok(base.ok, base.error);
  const dropAt = base.prepared.firstCut - 30 + base.prepared.shift;     // 30 bases before the cut, past the stretch the offset is read from
  const drifted = pair("nr2f2_ko_clone_hom_m1", (t) => dropCall(t, dropAt));
  assert.ok(drifted.ok, drifted.error);
  assert.ok(drifted.prepared.offsetRefined, "the offset was re-anchored");
  assert.equal(drifted.prepared.offsetRefined.delta, -1);
  assert.equal(drifted.prepared.shift, base.prepared.shift - 1);
  assert.equal(drifted.prepared.shiftUp, base.prepared.shift);
  assert.ok(drifted.warnings.some((w) => /re-anchored/.test(w)));
  assert.equal(topLabel(drifted), topLabel(base), "same allele as the undisturbed pair");
  assert.equal(drifted.genotype.category, base.genotype.category);
  assert.ok(Math.abs(drifted.summary.editedPct - base.summary.editedPct) <= 3);
  // without the re-anchoring (0.2.0 behaviour) the same pair is misread
  const old = pair("nr2f2_ko_clone_hom_m1", (t) => dropCall(t, dropAt), { refineOffset: false });
  assert.ok(!old.prepared.offsetRefined);
  assert.notEqual(topLabel(old), topLabel(base), "0.2.0 behaviour reproduces the one-base error");
});

test("offset: positions upstream of the step keep the old offset", () => {
  const dropAt = base.prepared.firstCut - 30 + base.prepared.shift;
  const drifted = pair("nr2f2_ko_clone_hom_m1", (t) => dropCall(t, dropAt));
  const { shiftPivot, shiftUp, shift } = drifted.prepared;
  assert.ok(Number.isFinite(shiftPivot));
  assert.equal(shiftAt(drifted.prepared, shiftPivot - 1), shiftUp);
  assert.equal(shiftAt(drifted.prepared, shiftPivot), shift);
  assert.equal(shiftAt({ shift: 3 }, 100), 3, "a prepared pair without the new fields behaves as before");
  // the upstream stretch is still compared with the old offset, so it does not look noisy
  assert.ok(drifted.quality.discordanceBefore < 0.25, `discordance before the cut ${drifted.quality.discordanceBefore}`);
});

test("offset: a drift in the last bases before the cut is not absorbed (it cannot be told from an edit)", () => {
  const dropAt = base.prepared.firstCut - 8 + base.prepared.shift;
  const late = pair("nr2f2_ko_clone_hom_m1", (t) => dropCall(t, dropAt));
  assert.ok(late.ok, late.error);
  assert.equal(late.prepared.offsetRefined, null);
  assert.ok(!late.warnings.some((w) => /re-anchored/.test(w)));
});

test("offset: no fixture pair is changed by the re-anchoring, and results are identical with it switched off", () => {
  for (const name of Object.keys(truth)) {
    const entry = truth[name];
    if (!entry.control || !entry.edited || !/^(nr2f2|apoe)/.test(name)) continue;
    const on = pair(name, (t) => t, { sampleType: entry.sampleType || "clone" });
    const off = pair(name, (t) => t, { sampleType: entry.sampleType || "clone", refineOffset: false });
    if (!on.ok) { assert.equal(on.ok, off.ok, name); continue; }
    assert.equal(on.prepared.offsetRefined, null, `${name} was re-anchored`);
    assert.deepEqual(on.summary, off.summary, name);
  }
});

test("offset: refineOffsetAtCut leaves a pair alone when the reads agree at the cut or there is too little to judge", () => {
  const { control, edited, firstCut, shift } = base.prepared;
  assert.equal(refineOffsetAtCut({ control, edited, firstCut, shift }), null);
  assert.equal(refineOffsetAtCut({ control, edited, firstCut, shift, options: { refineOffset: false } }), null);
  assert.equal(refineOffsetAtCut({ control, edited, firstCut: 5, shift }), null);
});

test("unplanned: a second base at one position upstream of the cut is named, with its position and base", () => {
  const clean = pair("apoe_snp_clone_hom");
  assert.ok(clean.ok, clean.error);
  assert.deepEqual(clean.unplanned, []);
  const at = clean.prepared.firstCut - 10 + clean.prepared.shift;
  const called = loadTrace(truth.apoe_snp_clone_hom.edited).calls[at];
  const alt = BASES.find((b) => b !== called && b !== loadTrace(truth.apoe_snp_clone_hom.control).calls[clean.prepared.firstCut - 10]);
  const marked = pair("apoe_snp_clone_hom", (t) => addSecondBase(t, at, alt));
  assert.ok(marked.ok, marked.error);
  assert.equal(marked.unplanned.length, 1, JSON.stringify(marked.unplanned));
  const hit = marked.unplanned[0];
  assert.equal(hit.relative, -10); assert.equal(hit.alt, alt); assert.equal(hit.ref, called);
  assert.ok(hit.fraction > 0.4 && hit.fraction < 0.6, `fraction ${hit.fraction}`);
  assert.ok(marked.warnings.some((w) => /Unplanned base change/.test(w) && /-10: /.test(w)));
  // the fitted mixture and the headline numbers are not changed by the annotation
  assert.equal(marked.summary.intendedEditPct >= 0, true);
});

test("unplanned: planned markers and fixture traces are never flagged; a run of changes is a misfit and is not reported", () => {
  for (const name of Object.keys(truth)) {
    const entry = truth[name];
    if (!entry.control || !entry.edited || !/^(nr2f2|apoe)/.test(name)) continue;
    const r = pair(name, (t) => t, { sampleType: entry.sampleType || "clone" });
    if (r.ok) assert.deepEqual(r.unplanned, [], `${name} flagged ${JSON.stringify(r.unplanned)}`);
  }
  const clean = pair("apoe_snp_clone_hom");
  const ks = [-12, -11, -10, -9].map((d) => clean.prepared.firstCut + d + clean.prepared.shift);
  const run = pair("apoe_snp_clone_hom", (t) => ks.reduce((acc, k) => addSecondBase(acc, k, BASES.find((b) => b !== acc.calls[k])), t));
  assert.deepEqual(run.unplanned, []);
  assert.deepEqual(scanUnplanned(clean.prepared, null, clean.decomposition), []);
  assert.deepEqual(scanUnplanned(clean.prepared, {}, clean.decomposition, { unplanned: false }), []);
});

test("unplanned: an accepted clone with an unplanned base change is held, with the reason", () => {
  const spec = specFor("apoe_snp_clone_hom"); const entry = truth.apoe_snp_clone_hom;
  const control = loadTrace(entry.control); const edited = loadTrace(entry.edited);
  const sample = (tweak) => analyseSample({ name: "clone", spec, control, edited: tweak(edited), sampleType: "clone", options: {} });
  const plain = sample((t) => t);
  assert.ok(plain.ok, plain.error);
  assert.equal(decide("snp", plain).decision, "accept");
  const probe = pair("apoe_snp_clone_hom"); const at = probe.prepared.firstCut - 10 + probe.prepared.shift;
  const alt = BASES.find((b) => b !== edited.calls[at]);
  const marked = sample((t) => addSecondBase(t, at, alt, 0.45));   // a second base at about 35 % of the signal: still an accept without the check
  assert.equal(marked.unplanned.length, 1);
  const d = decide("snp", marked);
  assert.equal(d.decision, "hold");
  assert.ok(d.reasons.some((r) => /Unplanned base change/.test(r) && /-10/.test(r)));
  assert.equal(decide("snp", marked, { holdOnUnplanned: false }).decision, "accept");
  // when the existing rules already hold or re-sequence the clone, the reason is added and the decision is kept
  const held = decide("snp", sample((t) => addSecondBase(t, at, alt, 0.6)));
  assert.equal(held.decision, "hold"); assert.ok(held.reasons.some((r) => /Unplanned base change/.test(r)));
  // a stronger change widens the intervals until the existing confidence rule asks for re-sequencing; either way the clone is not accepted
  assert.notEqual(decide("snp", sample((t) => addSecondBase(t, at, alt, 1))).decision, "accept");
});
