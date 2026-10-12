// A control with mixed signal (blurred peaks, e.g. after a homopolymer run) is explained, not reported as the wrong design.
import test from "node:test";
import assert from "node:assert/strict";
import { analysePair } from "../src/analyze.js";
import { describeMixedSignal, findHomopolymer } from "../src/diagnose.js";
import { peakHeights } from "../src/traceModel.js";
import { loadSpec, loadTrace } from "./helpers.js";

const BASES = ["A", "C", "G", "T"];
// Each peak becomes a weighted blend of its neighbours' peaks from base `from` on; the base calls follow the blended signal.
function blend(trace, from, weights = { "-2": 0.08, "-1": 0.25, 0: 0.34, 1: 0.25, 2: 0.08 }) {
  const heights = peakHeights(trace, 1);
  const n = trace.calls.length;
  const channels = Object.fromEntries(BASES.map((b) => [b, Int16Array.from(trace.channels[b])]));
  let calls = "";
  for (let i = 0; i < n; i += 1) {
    const mixed = BASES.map((_, b) => (i < from ? heights[i * 4 + b] : Object.entries(weights).reduce((sum, [d, w]) => sum + w * heights[Math.min(n - 1, Math.max(0, i + Number(d))) * 4 + b], 0)));
    const centre = trace.peakLocations[i];
    BASES.forEach((base, b) => { for (let s = centre - 1; s <= centre + 1; s += 1) if (s >= 0 && s < channels[base].length) channels[base][s] = Math.round(mixed[b]); });
    calls += BASES[mixed.indexOf(Math.max(...mixed))];
  }
  return { ...trace, channels, calls };
}

test("a control that turns mixed is called mixed, with the onset", () => {
  const spec = loadSpec("nr2f2_2xha_ssodn");
  const control = blend(loadTrace("nr2f2_control.ab1"), 150);
  const result = analysePair({ control, edited: loadTrace("nr2f2_ha_hdr_clone_het_edited.ab1"), spec, options: { sampleType: "clone" } });
  assert.equal(result.ok, false);
  assert.ok(result.mixedSignal, result.error);
  assert.match(result.error, /control trace itself is mixed/);
  assert.match(result.error, /not a wrong-design problem/);
  assert.match(result.error, /turns mixed from base \d+/);
});

test("an unrelated read is still reported as the wrong design", () => {
  const spec = loadSpec("apoe_r176c_snp");
  const result = analysePair({ control: loadTrace("nr2f2_control.ab1"), edited: loadTrace("nr2f2_ha_hdr_clone_het_edited.ab1"), spec });
  assert.equal(result.ok, false);
  assert.ok(!result.mixedSignal);
  assert.match(result.error, /does not match the design reference/);
});

test("the explanation names a homopolymer run just before the onset, and a clean read gets none", () => {
  const reference = "ACGT".repeat(30) + "A".repeat(20) + "GTCA".repeat(40);
  const n = 160; const purity = new Float64Array(n); const readToRef = new Int32Array(n);
  for (let k = 0; k < n; k += 1) { readToRef[k] = 100 + k; purity[k] = k < 50 ? 0.9 : 0.55; }
  const text = describeMixedSignal({ purity, readToRef, readStart: 0, readEnd: n, reference, identity: 0.74, alignedBases: n });
  assert.match(text, /20 bp A run/);
  assert.match(text, /agree at only 74%/);
  const clean = new Float64Array(n).fill(0.9);
  assert.equal(describeMixedSignal({ purity: clean, readToRef, readStart: 0, readEnd: n, reference, identity: 0.74, alignedBases: n }), null);
  assert.equal(findHomopolymer("ACGTAAAAAAAAAAAACGT", 0, 19).length, 12);
});
