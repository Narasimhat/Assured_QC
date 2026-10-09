import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { runJobs, parseBandsText } from "../src/runJobs.js";
import { loadSpec, loadTruth, simDir } from "./helpers.js";

const truth = loadTruth();
const file = (name) => { const bytes = readFileSync(path.join(simDir, name)); return { name, buffer: bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) }; };

test("bands text: named lines, a shared line, junk ignored", () => {
  const parsed = parseBandsText("clone1: 841, 631\nclone2: 631\n 700 ;800\n\nx: abc");
  assert.deepEqual(parsed.byName.get("clone1"), [841, 631]);
  assert.deepEqual(parsed.byName.get("clone2"), [631]);
  assert.deepEqual(parsed.all, [700, 800]);
  assert.deepEqual(parsed.byName.get("x"), []);
});

test("trace workflow with an exported design: every sample analysed, progress reported", () => {
  const seen = [];
  const out = runJobs({ workflow: "snp", designSpec: loadSpec("apoe_r176c_snp"), control: file("apoe_control.ab1"), samples: [file("apoe_snp_clone_het_edited.ab1"), file("apoe_snp_clone_partial_het_edited.ab1")], sampleType: "clone" }, (p) => seen.push(p.done));
  assert.ok(out.ok, out.error);
  assert.equal(out.source, "design");
  assert.deepEqual(out.results.map((r) => r.name), ["apoe_snp_clone_het_edited", "apoe_snp_clone_partial_het_edited"]);
  assert.deepEqual(out.results.map((r) => r.genotype.category), ["heterozygous_edit", "edit_plus_partial"]);
  assert.deepEqual(seen, [0, 1, 2]);
  assert.equal(out.design.markers, 5);
  JSON.stringify(out); // plain data: it crosses the worker boundary
});

test("manual entry: guides and donor typed in", () => {
  const entry = truth.apoe_snp_clone_het;
  const out = runJobs({ workflow: "snp", control: file(entry.control), samples: [file(entry.edited)], sampleType: "clone", manual: { guides: entry.guides.join(", "), donor: entry.donor, gene: "APOE" } });
  assert.ok(out.ok, out.error);
  assert.equal(out.source, "manual");
  assert.equal(out.results[0].genotype.category, "heterozygous_edit");
});

test("reporter workflow: junction reads from a design file, with out-out band sizes", () => {
  const out = runJobs({ workflow: "reporter", designSpec: loadSpec("nr2f2_sd40v5_block"), samples: [file("junction5_ok.ab1"), file("junction5_del_in_insert.ab1")], bandsText: "junction5_ok: 841\njunction5_del_in_insert: 841, 631" });
  assert.ok(out.ok, out.error);
  assert.deepEqual(out.results.map((r) => r.verdict.status), ["pass", "fail"]);
  assert.equal(out.results[0].bands.call, "knock-in only");
  assert.equal(out.results[1].bands.call, "knock-in and wild type");
  assert.match(runJobs({ workflow: "reporter", samples: [file("junction5_ok.ab1")] }).error, /need the exported design file/);
});

test("input problems are reported plainly", () => {
  assert.match(runJobs({ workflow: "nope", samples: [] }).error, /Unknown workflow/);
  assert.match(runJobs({ workflow: "knockout", samples: [] }).error, /at least one edited trace/);
  assert.match(runJobs({ workflow: "knockout", samples: [file("nr2f2_ko_pool_g1_edited.ab1")] }).error, /control/);
  assert.match(runJobs({ workflow: "snp", control: file("apoe_control.ab1"), samples: [file("apoe_snp_clone_het_edited.ab1")], manual: { guides: truth.apoe_snp_clone_het.guides.join(",") } }).error, /needs a donor/);
  const broken = runJobs({ workflow: "knockout", designSpec: loadSpec("nr2f2_2xha_ssodn"), control: file("nr2f2_control.ab1"), samples: [{ name: "notatrace.ab1", buffer: new ArrayBuffer(300) }, file("nr2f2_ko_pool_g1_edited.ab1")] });
  assert.ok(broken.ok);
  assert.equal(broken.results[0].ok, false);
  assert.match(broken.results[0].error, /Could not read notatrace\.ab1: This is not an ABIF/);
  assert.equal(broken.results[1].ok, true);
});
