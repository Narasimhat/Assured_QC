import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { analyseJunctionRead, classifyBands } from "../src/junction.js";
import { loadSpec, loadTrace, simDir } from "./helpers.js";

const spec = loadSpec("nr2f2_sd40v5_block");
const truth = JSON.parse(readFileSync(path.join(simDir, "junction_truth.json"), "utf8"));
const read = (name) => analyseJunctionRead({ trace: loadTrace(truth[name].file), spec });

test("a correct 5' junction read passes and shows the insert, both junction sides of the donor and the genomic flank", () => {
  const result = read("junction5_ok");
  assert.ok(result.ok, result.error);
  assert.equal(result.verdict.status, "pass", result.verdict.details.join(" | "));
  assert.equal(result.orientation, "+");
  assert.ok(result.identity > 0.99);
  assert.equal(result.spans.fivePrimeJunction, true);
  assert.equal(result.spans.threePrimeJunction, true, "this long read crosses the whole cassette");
  assert.ok(result.spans.genomicFlankBases >= 100);
  assert.equal(result.events.filter((e) => e.type !== "mismatch").length, 0);
  const covered = result.markers.filter((m) => m.covered);
  assert.ok(covered.length >= 1 && covered.every((m) => m.present), "the blocking change in the 5' arm is present");
});

test("a read from the reverse primer across the 3' junction is oriented and passes", () => {
  const result = read("junction3_ok");
  assert.ok(result.ok, result.error);
  assert.equal(result.orientation, "-");
  assert.equal(result.verdict.status, "pass", result.verdict.details.join(" | "));
  assert.equal(result.spans.threePrimeJunction, true);
});

test("a one-base deletion inside the insert fails and is placed in the insert", () => {
  const result = read("junction5_del_in_insert");
  assert.ok(result.ok, result.error);
  assert.equal(result.verdict.status, "fail");
  const deletion = result.events.find((e) => e.type === "deletion");
  assert.ok(deletion, "deletion reported");
  assert.equal(deletion.region, "insert");
  assert.equal(deletion.length, 1);
  const meta = truth._meta;
  assert.ok(Math.abs(deletion.alleleIndex - (meta.insert_start_in_allele + 100)) <= 2, `at ${deletion.alleleIndex}`);
});

test("a substitution in the 5' arm is a warning with its location", () => {
  const result = read("junction5_arm_snp");
  assert.ok(result.ok, result.error);
  assert.equal(result.verdict.status, "warn", result.verdict.details.join(" | "));
  const mismatch = result.events.find((e) => e.type === "mismatch" && !e.expectedChange);
  assert.equal(mismatch.region, "5' arm");
  assert.equal(result.events.filter((e) => e.type !== "mismatch").length, 0);
});

test("a missing blocking change is reported as not copied", () => {
  const result = read("junction5_blocking_missing");
  assert.ok(result.ok, result.error);
  assert.equal(result.verdict.status, "warn");
  assert.ok(result.verdict.details.some((text) => /absent: the read shows the original base/.test(text)));
  const marker = result.markers.find((m) => m.covered && !m.present);
  assert.ok(marker && marker.reverted);
});

test("a wild-type read and an unrelated read are refused with the right reason", () => {
  const wild = analyseJunctionRead({ trace: loadTrace("nr2f2_control.ab1"), spec });
  assert.equal(wild.ok, false);
  assert.match(wild.error, /matches the wild-type sequence at least as well/);
  assert.equal(wild.looksWildType, true);
  const other = analyseJunctionRead({ trace: loadTrace("apoe_control.ab1"), spec });
  assert.equal(other.ok, false);
  assert.match(other.error, /does not match the expected knock-in allele/);
  assert.equal(analyseJunctionRead({ trace: loadTrace("nr2f2_control.ab1"), spec: { ...spec, donors: [] } }).ok, false);
});

test("out-out PCR band sizes: wild type, knock-in, both, unexpected", () => {
  const amplicon = spec.amplicons[0];
  const ki = Object.values(amplicon.editedBp)[0];
  assert.deepEqual([amplicon.wtBp, ki], [631, 841]);
  assert.equal(classifyBands([841], amplicon).call, "knock-in only");
  assert.equal(classifyBands([631], amplicon).call, "wild type only");
  assert.equal(classifyBands([850, 640], amplicon).call, "knock-in and wild type");
  assert.equal(classifyBands([841, 500], amplicon).call, "unexpected sizes");
  assert.equal(classifyBands([], amplicon).call, "no band");
  assert.equal(classifyBands([841], amplicon).bands[0].matches, "knock-in");
  assert.equal(classifyBands([1400], amplicon).call, "unexpected sizes");
});
