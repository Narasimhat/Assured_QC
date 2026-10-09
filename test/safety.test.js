// Misuse and failure cases: what the engine must refuse, warn about, or leave at zero.
import test from "node:test";
import assert from "node:assert/strict";
import { analysePair, NOISY_DISCORDANCE } from "../src/analyze.js";
import { runJobs } from "../src/runJobs.js";
import { readFileSync } from "node:fs";
import path from "node:path";
import { loadSpec, loadTrace, simDir } from "./helpers.js";

const spec = loadSpec("nr2f2_2xha_ssodn");
const control = loadTrace("nr2f2_control.ab1");
const het = loadTrace("nr2f2_ha_hdr_clone_het_edited.ab1");
const analyse = (edited, ctrl = control) => analysePair({ control: ctrl, edited, spec, options: { sampleType: "clone" } });

const truncate = (trace, samples) => {
  const keep = trace.peakLocations.filter((p) => p <= samples).length;
  return { ...trace, samples, channels: Object.fromEntries(Object.entries(trace.channels).map(([c, v]) => [c, v.slice(0, samples)])), calls: trace.calls.slice(0, keep), peakLocations: trace.peakLocations.slice(0, keep), quality: trace.quality.slice(0, keep) };
};
// Deterministic Gaussian noise added to every channel.
const noisy = (trace, sd) => {
  let seed = 7; const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
  return { ...trace, channels: Object.fromEntries(Object.entries(trace.channels).map(([c, v]) => [c, Int16Array.from(v, (x) => Math.max(0, Math.round(x + gauss() * sd)))])) };
};

test("the control analysed as its own edited sample shows no editing", () => {
  const result = analyse(control);
  assert.ok(result.ok, result.error);
  assert.equal(result.summary.editedPct, 0);
  assert.equal(result.genotype.category, "homozygous_wt");
});

test("control and edited traces swapped: refused, not reported as an edit", () => {
  const result = analyse(control, het);
  assert.equal(result.ok, false);
  assert.match(result.error, /control read does not match the design reference/);
});

test("a read that stops before or just after the cut is refused with a reason", () => {
  const early = analyse(truncate(het, Math.round(het.samples * 0.15)));
  assert.equal(early.ok, false);
  assert.match(early.error, /does not align to the control upstream of the cut/);
  const short = analyse(truncate(het, Math.round(het.samples * 0.4)));
  assert.equal(short.ok, false);
  assert.match(short.error, /too short/);
});

test("moderate noise is accepted silently; heavy noise raises a warning", () => {
  const moderate = analyse(noisy(het, 50));
  assert.ok(moderate.ok, moderate.error);
  assert.ok(!moderate.warnings.some((w) => /already differ before the cut/.test(w)), moderate.warnings.join(" | "));
  assert.ok(Math.abs(moderate.summary.intendedEditPct - 50) < 6, `intended ${moderate.summary.intendedEditPct}`);
  const heavy = analyse(noisy(het, 250));
  assert.ok(heavy.ok, heavy.error);
  assert.ok(heavy.quality.discordanceBefore > NOISY_DISCORDANCE, `discordance ${heavy.quality.discordanceBefore}`);
  assert.ok(heavy.warnings.some((w) => /already differ before the cut/.test(w)));
});

test("a damaged file in a batch is reported by name and the other samples are still analysed", () => {
  const fixture = (name) => { const bytes = readFileSync(path.join(simDir, name)); return { name, buffer: bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) }; };
  const damaged = { name: "bad.ab1", buffer: Uint8Array.from([1, 2, 3, 4, 5]).buffer };
  const out = runJobs({ workflow: "knockout", designSpec: spec, control: fixture("nr2f2_control.ab1"), samples: [damaged, fixture("nr2f2_ko_pool_g1_edited.ab1")], sampleType: "pool" });
  assert.ok(out.ok, out.error);
  assert.equal(out.results[0].ok, false);
  assert.match(out.results[0].error, /Could not read bad\.ab1/);
  assert.equal(out.results[1].ok, true);
  assert.ok(out.results[1].summary.editedPct > 60);
});
