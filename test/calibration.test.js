import test from "node:test";
import assert from "node:assert/strict";
import { wilson, parseObservations, calibrate, calibrationMarkdown } from "../src/calibration.js";

const result = (name, status, cls, fraction = 1) => ({ name, summary: { functionalKoPct: 90 }, functional: { call: { status }, rows: [{ label: "allele", class: cls, fraction, lof: status.startsWith("predicted") }, { label: "wt", class: "wt", fraction: 0, lof: false }] } });

test("Wilson interval matches known values", () => {
  const w = wilson(8, 10); assert.ok(Math.abs(w.p - 0.8) < 1e-12); assert.ok(Math.abs(w.lo - 0.4902) < 5e-4 && Math.abs(w.hi - 0.9433) < 5e-4, `${w.lo} ${w.hi}`);
  assert.deepEqual(wilson(0, 0), { p: null, lo: null, hi: null }); assert.equal(wilson(10, 10).hi, 1);
});

test("observations are read from comma or tab separated text and bad input is reported", () => {
  const csv = parseObservations("sample,protein,rna\nA1,absent,decayed\n# comment\nA2,\"present\",stable\n"); assert.equal(csv.errors.length, 0); assert.equal(csv.rows.length, 2); assert.equal(csv.rows[1].protein, "present");
  assert.equal(parseObservations("sample\tprotein\nA1\tabsent\n").rows[0].protein, "absent");
  assert.ok(parseObservations("clone,western\nA1,absent\n").errors.length >= 1);
});

test("calls are scored against protein and the failing allele classes are found", () => {
  const results = [result("c1", "predicted_knockout", "lof_nmd"), result("c2", "predicted_knockout", "lof_nmd"), result("c3", "predicted_knockout_with_caveats", "lof_truncation"), result("c4", "predicted_knockout", "lof_nmd"), result("c5", "uncertain", "inframe"), result("c6", "not_knockout", "wt"), result("c7", "uncertain", "truncation_uncertain")];
  const obs = [["c1", "absent", "decayed"], ["c2", "absent", "decayed"], ["c3", "present", ""], ["c4", "reduced", "stable"], ["c5", "present", ""], ["c6", "present", ""], ["c7", "absent", ""], ["zz", "absent", ""]].map(([sample, protein, rna]) => ({ sample, protein, rna }));
  const cal = calibrate(results, obs);
  assert.equal(cal.n, 7); assert.equal(cal.problems.length, 1); assert.match(cal.problems[0], /zz/);
  const k = cal.knockoutCall; assert.equal(k.called, 4); assert.equal(k.confirmed, 2); assert.equal(k.notConfirmed, 2); assert.equal(k.notCalled, 3); assert.equal(k.missed, 1);
  assert.ok(Math.abs(k.sensitivity.p - 2 / 3) < 1e-12); assert.ok(Math.abs(k.specificity.p - 2 / 4) < 1e-12, "tn=2 of 4 without the knockout call and without protein loss");
  const nmd = cal.byClass.find((c) => c.class === "lof_nmd"); assert.equal(nmd.n, 3); assert.equal(nmd.absent, 2);
  assert.equal(cal.nmd.withRna, 3); assert.equal(cal.nmd.decayed, 2);
  assert.match(cal.note, /Only 7 clones/);
  const md = calibrationMarkdown(cal); assert.match(md, /predicted_knockout \| 3 \| 2/); assert.match(md, /Not matched/);
});

test("a protein word that is not recognised is reported, not guessed", () => {
  const cal = calibrate([result("c1", "predicted_knockout", "lof_nmd")], [{ sample: "c1", protein: "faint smear" }]);
  assert.equal(cal.n, 0); assert.match(cal.problems[0], /not absent, reduced, present or truncated/);
});
