// Expected alleles are built from the design alone, then checked against allele sequences that an
// independent simulator (Python, scripts/make_qc_fixtures.py) used to build the synthetic traces.
import test from "node:test";
import assert from "node:assert/strict";
import { buildExpectedAlleles } from "../src/alleleLibrary.js";
import { reverseSpec } from "../src/designSpec.js";
import { loadSpec, loadTruth } from "./helpers.js";

const truth = loadTruth();
const byKind = (library, kind) => library.alleles.filter((allele) => allele.kind === kind);

test("wild type and unique sequences", () => {
  for (const name of ["nr2f2_sd40v5_block", "apoe_r176c_snp", "tagme_ko_deletion"]) {
    const spec = loadSpec(name);
    const library = buildExpectedAlleles(spec);
    assert.equal(library.alleles[0].kind, "wt");
    assert.equal(library.alleles[0].sequence, spec.reference.sequence);
    assert.equal(new Set(library.alleles.map((a) => a.sequence)).size, library.alleles.length, `${name}: duplicate sequences`);
    assert.equal(new Set(library.alleles.map((a) => a.id)).size, library.alleles.length, `${name}: duplicate ids`);
  }
});

test("terminal tag block: three features (two blocking changes and the insert) give six conversion runs", () => {
  const spec = loadSpec("nr2f2_sd40v5_block");
  const library = buildExpectedAlleles(spec);
  const donor = spec.donors[0];
  assert.equal(library.donorAlleles.length, 6);
  const complete = byKind(library, "edit");
  assert.equal(complete.length, 1);
  assert.equal(complete[0].sequence, spec.reference.sequence.slice(0, donor.refStart) + donor.sequence + spec.reference.sequence.slice(donor.refEnd));
  assert.equal(complete[0].size, donor.sequence.length - (donor.refEnd - donor.refStart));
  assert.equal(complete[0].size, 210);
  // the insert without any blocking change is one of the partial runs
  assert.ok(library.donorAlleles.some((a) => a.carries.length === 1 && a.carries[0] === "insert"));
});

test("SNP with silent blocking changes: every contiguous run is an allele, and the silent-only run lacks the SNP", () => {
  const spec = loadSpec("apoe_r176c_snp");
  const library = buildExpectedAlleles(spec);
  assert.equal(library.donorAlleles.length, 2 * 18, "15 contiguous runs plus 3 all-but-one sets per donor");
  const donor = spec.donors[0];
  const intended = spec.markers.find((marker) => marker.role === "intended");
  const silentOnly = library.donorAlleles.filter((allele) => allele.donor === donor.name && allele.kind === "edit_partial" && !allele.intendedPresent);
  assert.ok(silentOnly.length >= 4, "runs that exclude the SNP exist");
  const all = library.donorAlleles.find((allele) => allele.donor === donor.name && allele.kind === "edit");
  const withoutSnp = all.sequence.slice(0, intended.pos) + intended.ref + all.sequence.slice(intended.pos + 1);
  assert.ok(silentOnly.some((allele) => allele.sequence === withoutSnp), "all blocking changes but not the SNP is present");
  assert.ok(byKind(library, "edit").every((allele) => allele.size === 0));
});

test("indels cover every deletion size at each cut, and short insertions", () => {
  const spec = loadSpec("nr2f2_sd40v5_block");
  const library = buildExpectedAlleles(spec, { maxDeletion: 12 });
  spec.guides.forEach((guide, index) => {
    const sizes = new Set(byKind(library, "indel").filter((a) => a.guide === index + 1 && a.size < 0).map((a) => -a.size));
    assert.deepEqual([...sizes].sort((a, b) => a - b), Array.from({ length: 12 }, (_, k) => k + 1));
    const plus1 = byKind(library, "indel").filter((a) => a.guide === index + 1 && a.size === 1);
    assert.equal(plus1.length, 4);
    plus1.forEach((a) => assert.equal(a.sequence.length, spec.reference.sequence.length + 1));
  });
  assert.ok(byKind(library, "indel").every((a) => a.sequence.length === spec.reference.sequence.length + a.size));
});

test("two guides: the deletion between the cuts is an allele", () => {
  const spec = loadSpec("nr2f2_sd40v5_block");
  const library = buildExpectedAlleles(spec);
  const [g1, g2] = spec.guides;
  const exact = spec.reference.sequence.slice(0, g1.cut) + spec.reference.sequence.slice(g2.cut);
  const hit = byKind(library, "deletion_between_cuts").find((a) => a.sequence === exact);
  assert.ok(hit);
  assert.equal(hit.size, -(g2.cut - g1.cut));
  assert.equal(g2.cut - g1.cut, 29);
});

test("every allele used to build the synthetic traces is in the library (independent code paths)", () => {
  const specFor = { nr2f2: loadSpec("nr2f2_2xha_ssodn"), apoe: loadSpec("apoe_r176c_snp") };
  let checked = 0;
  for (const [name, entry] of Object.entries(truth)) {
    let spec = specFor[name.startsWith("nr2f2") ? "nr2f2" : "apoe"];
    // a read from the reverse primer is built on the reverse complement of the window
    if (spec.reference.sequence.indexOf(entry.reference) < 0) spec = reverseSpec(spec);
    const library = buildExpectedAlleles(spec);
    const start = spec.reference.sequence.indexOf(entry.reference);
    assert.ok(start >= 0, `${name}: simulated read lies inside the design window`);
    const span = 700;
    for (const allele of entry.alleles) {
      const hit = library.alleles.find((candidate) => candidate.sequence.slice(start, start + span) === allele.sequence.slice(0, span));
      assert.ok(hit, `${name}: allele ${allele.label} missing from the library`);
      checked += 1;
    }
  }
  assert.equal(checked, Object.values(truth).reduce((sum, entry) => sum + entry.alleles.length, 0));
});

test("deletions are placed so that they remove the cut: from starting at it to ending at it, and no further", () => {
  const spec = loadSpec("nr2f2_sd40v5_block");
  const reference = spec.reference.sequence;
  const { cut } = spec.guides[0];
  const library = buildExpectedAlleles(spec, { maxDeletion: 5 });
  const have = new Set(library.alleles.map((allele) => allele.sequence));
  const size = 5;
  for (let offset = 0; offset <= size; offset += 1) {
    const start = cut - offset;
    assert.ok(have.has(reference.slice(0, start) + reference.slice(start + size)), `deletion of ${size} starting ${offset} bp before the cut`);
  }
  const valid = new Set(Array.from({ length: size + 1 }, (_, offset) => reference.slice(0, cut - offset) + reference.slice(cut - offset + size)));
  const tooFarRight = reference.slice(0, cut + 1) + reference.slice(cut + 1 + size);
  const tooFarLeft = reference.slice(0, cut - size - 1) + reference.slice(cut - 1);
  [tooFarRight, tooFarLeft].forEach((sequence) => { if (!valid.has(sequence)) assert.ok(!have.has(sequence) || spec.guides.length > 1, "a deletion that misses the cut is not in the library"); });
});
