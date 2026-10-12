import test from "node:test";
import assert from "node:assert/strict";
import { analyseReads, analysePair } from "../src/analyze.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const truth = loadTruth(); const spec = loadSpec("nr2f2_2xha_ssodn");
const read = (name, label) => ({ name: label, control: loadTrace(truth[name].control), edited: loadTrace(truth[name].edited) });
const wtShare = (name) => truth[name].alleles.filter((a) => a.label === "wt").reduce((s, a) => s + a.fraction, 0) * 100;

test("forward and reverse reads of the same pool are fitted together and agree", () => {
  const joint = analyseReads({ reads: [read("nr2f2_ko_pool_g1", "F"), read("nr2f2_rev_ko_pool_g1", "R")], spec, options: { sampleType: "pool" } });
  assert.ok(joint.ok, joint.error); assert.equal(joint.readsUsed, 2);
  const edited = 100 - wtShare("nr2f2_ko_pool_g1");
  assert.ok(Math.abs(joint.summary.editedPct - edited) < 6, `joint ${joint.summary.editedPct} vs truth ${edited}`);
  assert.equal(joint.reads.length, 2);
  assert.ok(joint.agreement.agree, JSON.stringify(joint.reads));
  const [lo, hi] = joint.intervals.editedPct; assert.ok(lo <= edited && edited <= hi, `${lo}-${hi} vs ${edited}`);
});

test("one read that cannot be used is reported and the other still gives a result", () => {
  const bad = { name: "R", control: loadTrace(truth.nr2f2_rev_ko_pool_g1.control), edited: loadTrace("apoe_control.ab1") };
  const result = analyseReads({ reads: [read("nr2f2_ko_pool_g1", "F"), bad], spec, options: { sampleType: "pool" } });
  assert.ok(result.ok, result.error); assert.equal(result.readsUsed, 1);
  assert.ok(result.warnings.some((w) => /Read R was not used/.test(w)), result.warnings.join(" | "));
});

test("reads of different samples disagree and the disagreement is flagged, not averaged away", () => {
  const mismatched = analyseReads({ reads: [read("nr2f2_ko_clone_hom_m1", "F"), { ...read("nr2f2_rev_ko_pool_g1", "R") }], spec, options: { sampleType: "pool" } });
  assert.ok(mismatched.ok, mismatched.error);
  assert.equal(mismatched.agreement.agree, false, JSON.stringify(mismatched.reads));
  assert.ok(mismatched.warnings.some((w) => /disagree/.test(w)));
});

test("a single read through analyseReads equals analysePair", () => {
  const one = analyseReads({ reads: [read("nr2f2_ko_clone_het_wt_m1", "F")], spec, options: { sampleType: "clone" } });
  const pair = analysePair({ control: loadTrace(truth.nr2f2_ko_clone_het_wt_m1.control), edited: loadTrace(truth.nr2f2_ko_clone_het_wt_m1.edited), spec, options: { sampleType: "clone" } });
  assert.equal(one.summary.editedPct, pair.summary.editedPct);
});
