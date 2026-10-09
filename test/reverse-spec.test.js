import test from "node:test";
import assert from "node:assert/strict";
import { readdirSync } from "node:fs";
import { reverseSpec, validateDesignSpec } from "../src/designSpec.js";
import { buildExpectedAlleles } from "../src/alleleLibrary.js";
import { reverseComplement } from "../src/seq.js";
import { designDir, loadSpec } from "./helpers.js";

const names = readdirSync(designDir).map((file) => file.replace(".design.json", ""));

test("the reversed design is still a valid design", () => {
  for (const name of names) {
    const reversed = reverseSpec(loadSpec(name));
    const { errors } = validateDesignSpec({ ...reversed, amplicons: [] });
    assert.deepEqual(errors, [], `${name}: ${errors.join("; ")}`);
  }
});

test("reversing twice returns the original design", () => {
  for (const name of names) {
    const spec = loadSpec(name);
    const twice = reverseSpec(reverseSpec(spec));
    assert.equal(twice.reference.sequence, spec.reference.sequence);
    assert.deepEqual(twice.guides.map((g) => [g.strand, g.cut, g.protospacerStart, g.protospacerEnd]), spec.guides.map((g) => [g.strand, g.cut, g.protospacerStart, g.protospacerEnd]));
    assert.deepEqual(twice.markers.map((m) => [m.pos, m.ref, m.alt]), spec.markers.map((m) => [m.pos, m.ref, m.alt]));
    assert.deepEqual(twice.donors.map((d) => [d.refStart, d.refEnd, d.sequence, d.arm5, d.arm3, d.insertStart]), spec.donors.map((d) => [d.refStart, d.refEnd, d.sequence, d.arm5, d.arm3, d.insertStart]));
  }
});

test("alleles built on the reversed design are the reverse complements of the original alleles", () => {
  for (const name of ["nr2f2_sd40v5_block", "apoe_r176c_snp", "tagme_it_spot", "tagme_nt_egfp"]) {
    const spec = loadSpec(name);
    const forward = new Set(buildExpectedAlleles(spec).alleles.map((allele) => reverseComplement(allele.sequence)));
    const reversed = buildExpectedAlleles(reverseSpec(spec)).alleles.map((allele) => allele.sequence);
    assert.equal(reversed.length, forward.size, name);
    assert.ok(reversed.every((sequence) => forward.has(sequence)), `${name}: a reversed allele has no counterpart`);
  }
});
