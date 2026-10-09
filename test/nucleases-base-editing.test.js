import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { locateGuide, NUCLEASES } from "../src/nucleases.js";
import { baseEditMarkers } from "../src/baseEdit.js";
import { discoverCut } from "../src/cutDiscovery.js";
import { planDesign } from "../src/batch.js";
import { runJobs } from "../src/runJobs.js";
import { loadTruth, loadTrace } from "./helpers.js";

const here = path.dirname(fileURLToPath(import.meta.url));
const extended = JSON.parse(readFileSync(path.join(here, "..", "fixtures", "sim", "truth_extended.json"), "utf8"));
const truth = loadTruth();
const rc = (s) => [...s].reverse().map((b) => ({ A: "T", C: "G", G: "C", T: "A" }[b])).join("");

test("SpCas9 cuts 3 bases from the PAM on either strand; SaCas9 and Cas12a use their own PAM and cut", () => {
  const spacer = "GATTACAGATTACAGATTAC"; const plus = `CCCCCCCCCC${spacer}TGGAAAAAAAAAA`;
  const sp = locateGuide(plus, spacer, "SpCas9"); assert.equal(sp.strand, "+"); assert.equal(sp.pam, "TGG"); assert.ok(sp.pamOk); assert.equal(sp.cut, 10 + 20 - 3);
  const minus = rc(plus); const sm = locateGuide(minus, spacer, "SpCas9"); assert.equal(sm.strand, "-"); assert.equal(sm.pam, "TGG"); assert.ok(sm.pamOk); assert.equal(sm.cut, minus.indexOf(rc(spacer)) + 3);
  const sa = "GATTACAGATTACAGATTACAA"; const saRef = `CCCCCCCCCC${sa}TTGAGTAAAAAAA`; const s2 = locateGuide(saRef, sa, "SaCas9"); assert.equal(s2.pam, "TTGAGT"); assert.ok(s2.pamOk);
  const cas12 = `CCCCCCCCCCTTTA${spacer}AAAAAAAAAA`; const c12 = locateGuide(cas12, spacer, "Cas12a"); assert.equal(c12.pam, "TTTA"); assert.ok(c12.pamOk); assert.equal(c12.cut, 14 + 20);
  const wrong = locateGuide(plus, spacer, "Cas12a"); assert.equal(wrong.pamOk, false);
  assert.ok(Object.keys(NUCLEASES).length >= 5);
});

test("base-editing markers: ABE and CBE positions in the window, on either strand", () => {
  const spacer = "ACACACACACACACACACAC"; const reference = `GGGGGGGGGG${spacer}CGGTTTTTTTTT`;
  const abe = baseEditMarkers({ reference, guide: spacer, editor: "ABE", window: [4, 8], target: 6 });
  assert.deepEqual(abe.markers.map((m) => m.protospacerPosition), [5, 7]); assert.ok(abe.markers.every((m) => m.ref === "A" && m.alt === "G"));
  const cbe = baseEditMarkers({ reference, guide: spacer, editor: "CBE", window: [4, 8], target: 4 });
  assert.deepEqual(cbe.markers.map((m) => m.protospacerPosition), [4, 6, 8]); assert.equal(cbe.markers[0].role, "intended");
  const minusRef = rc(reference); const abeMinus = baseEditMarkers({ reference: minusRef, guide: spacer, editor: "ABE", window: [4, 8] });
  assert.deepEqual(abeMinus.markers.map((m) => m.protospacerPosition), [5, 7]); assert.ok(abeMinus.markers.every((m) => m.ref === "T" && m.alt === "C"));
});

test("cut discovery finds where an edited trace departs from its control", () => {
  const entry = truth.nr2f2_ko_clone_hom_m1; const found = discoverCut({ control: loadTrace(entry.control), edited: loadTrace(entry.edited) });
  assert.ok(found.ok, found.error); assert.ok(Math.abs(found.readIndex - entry.cuts.g1) <= 4, `found ${found.readIndex}, cut ${entry.cuts.g1}`);
  const none = discoverCut({ control: loadTrace(extended.ext_nr2f2_wt_replicate.control), edited: loadTrace(extended.ext_nr2f2_wt_replicate.edited) });
  assert.equal(none.ok, false);
});

test("with no guide on record the cut is inferred from the edited traces, and the analysis runs on it", () => {
  const entry = truth.nr2f2_ko_clone_het_wt_m1; const control = loadTrace(entry.control);
  const plan = planDesign({ control, guides: "", donor: "", gene: "x", editedForCut: [loadTrace(entry.edited)] });
  assert.ok(!plan.error, plan.error); assert.equal(plan.source, "inferred cut"); assert.ok(Math.abs(plan.spec.guides[0].cut - entry.cuts.g1) <= 4);
  assert.ok(plan.warnings.some((w) => /inferred/.test(w)));
  const none = planDesign({ control, guides: "", donor: "", gene: "x", editedForCut: [] }); assert.ok(none.error);
});

test("base editing: the per-position conversion is read from the traces", () => {
  const entry = extended.ext_nr2f2_abe_p8_40; const file = (name) => { const buffer = readFileSync(path.join(here, "..", "fixtures", "sim", name)); return { name, buffer: buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength) }; };
  const result = runJobs({ workflow: "base_edit", manual: { guides: entry.guides.join(","), gene: "nr2f2" }, baseEditor: { editor: "ABE", window: [4, 8], target: 8 }, sampleType: "pool", control: file(entry.control), samples: [file(entry.edited)] });
  assert.ok(result.ok, result.error); const sample = result.results[0]; assert.ok(sample.ok, sample.error);
  const row = sample.baseEdits.positions.find((p) => p.protospacerPosition === 8); assert.ok(row && row.covered);
  // the idealised fixture has random per-peak heights (12%), which the context correction cannot know; real controls have systematic ones
  assert.ok(Math.abs(row.conversionPct - 40) < 11, `conversion ${row.conversionPct} (uncorrected ${row.naivePct})`); assert.ok(Math.abs(row.naivePct - 40) < 6); assert.ok(row.intervalPct[0] <= 40 && row.intervalPct[1] >= 40);
  assert.equal(row.role, "intended");
});
