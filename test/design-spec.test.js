// The design hand-off must agree with itself before anything is analysed.
import test from "node:test";
import assert from "node:assert/strict";
import { readdirSync } from "node:fs";
import { validateDesignSpec } from "../src/designSpec.js";
import { designDir, loadSpec } from "./helpers.js";
const names = readdirSync(designDir).map((file) => file.replace(".design.json", ""));

test("every exported design validates", () => {
  assert.equal(names.length, 7);
  for (const name of names) {
    const result = validateDesignSpec(loadSpec(name));
    assert.deepEqual(result.errors, [], `${name}: ${result.errors.join("; ")}`);
  }
});

test("guide cut, protospacer and strand agree with the reference", () => {
  const spec = loadSpec("nr2f2_2xha_ssodn_oppstrand");
  assert.deepEqual(spec.guides.map((guide) => guide.strand), ["+", "-"]);
  assert.ok(validateDesignSpec(spec).ok);
});

test("a corrupted hand-off is rejected, one defect at a time", () => {
  const mutate = (name, edit) => { const spec = loadSpec(name); edit(spec); return validateDesignSpec(spec); };
  const base = "apoe_r176c_snp";
  assert.ok(mutate(base, (s) => { s.markers[0].alt = s.markers[0].ref; }).errors.some((e) => e.includes("alternate base equals")));
  assert.ok(mutate(base, (s) => { s.markers[1].ref = s.markers[1].ref === "A" ? "C" : "A"; }).errors.some((e) => e.includes("reference base is")));
  assert.ok(mutate(base, (s) => { s.guides[0].cut += 1; }).errors.some((e) => e.includes("not 3 bp from the PAM")));
  assert.ok(mutate(base, (s) => { s.guides[0].spacer = s.guides[0].spacer.replace(/^./, (b) => (b === "A" ? "C" : "A")); }).errors.some((e) => e.includes("does not match the reference")));
  assert.ok(mutate(base, (s) => { s.donors[0].sequence = s.donors[0].sequence.slice(1); }).errors.some((e) => e.includes("does not equal the reference span")));
  assert.ok(mutate(base, (s) => { s.markers.pop(); }).errors.some((e) => e.includes("no marker declares it")));
  assert.ok(mutate(base, (s) => { s.primers[0].sequence = "ACGTACGTACGTACGTAC"; }).errors.some((e) => e.includes("does not match the reference")));
  assert.ok(mutate(base, (s) => { s.amplicons[0].wtBp += 1; }).errors.some((e) => e.includes("stated wild-type size")));
  assert.ok(mutate(base, (s) => { s.schema = "other/9"; }).errors.some((e) => e.includes("Unsupported schema")));
  assert.ok(mutate(base, (s) => { s.reference.sequence = s.reference.sequence.replace(/A/, "R"); }).errors.some((e) => e.includes("A, C, G and T only")));
  assert.equal(validateDesignSpec(null).ok, false);
});

test("amplicon sizes carry the net insertion of each donor", () => {
  const spec = loadSpec("nr2f2_sd40v5_block");
  const donor = spec.donors[0];
  const amplicon = spec.amplicons[0];
  assert.equal(amplicon.editedBp[spec.donors[0].name === donor.name ? Object.keys(amplicon.editedBp)[0] : ""], amplicon.wtBp + donor.sequence.length - (donor.refEnd - donor.refStart));
  assert.deepEqual([amplicon.wtBp, Object.values(amplicon.editedBp)[0]], [631, 841]);
});
