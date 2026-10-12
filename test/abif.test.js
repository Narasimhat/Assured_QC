import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { parseAbif } from "../src/abif.js";
import { simDir, loadTruth } from "./helpers.js";

const load = (name) => parseAbif(readFileSync(path.join(simDir, name)), name);

test("parses channels, base calls, peak locations and quality from a trace", () => {
  const trace = load("nr2f2_control.ab1");
  assert.equal(trace.order, "GATC");
  assert.equal(trace.calls.length, 760);
  assert.equal(trace.peakLocations.length, 760);
  assert.equal(trace.quality.length, 760);
  assert.ok(trace.hasQuality);
  assert.ok(["A", "C", "G", "T"].every((base) => trace.channels[base].length === trace.samples));
  assert.ok(trace.peakLocations.every((p, i) => i === 0 || p > trace.peakLocations[i - 1]), "peak locations increase");
});

test("the called base is the tallest channel at its peak (control trace)", () => {
  const trace = load("apoe_control.ab1");
  let agree = 0;
  trace.peakLocations.forEach((p, i) => {
    const tallest = ["A", "C", "G", "T"].reduce((best, base) => (trace.channels[base][p] > trace.channels[best][p] ? base : best), "A");
    if (tallest === trace.calls[i]) agree += 1;
  });
  assert.ok(agree / trace.calls.length > 0.99, `${agree}/${trace.calls.length}`);
});

test("the control's calls are the design reference they were generated from", () => {
  const truth = loadTruth();
  const trace = load("nr2f2_control.ab1");
  assert.equal(trace.calls, truth.nr2f2_ko_pool_g1.reference);
});

test("damaged and non-ab1 inputs are refused with a clear message", () => {
  assert.throws(() => parseAbif(new Uint8Array(10)), /too small/);
  assert.throws(() => parseAbif(new Uint8Array(200)), /not an ABIF/);
  const bytes = new Uint8Array(readFileSync(path.join(simDir, "nr2f2_control.ab1")));
  const broken = bytes.slice(); broken[26] = 0x7f; // directory offset beyond the file
  assert.throws(() => parseAbif(broken), /directory is damaged|outside the file/);
  const noChannels = bytes.slice(); // rename DATA9 entry so the channel is missing
  const text = new TextDecoder("latin1").decode(noChannels);
  const at = text.lastIndexOf("DATA\u0000\u0000\u0000\u0009");
  assert.ok(at > 0);
  noChannels[at + 7] = 0x63;
  assert.throws(() => parseAbif(noChannels), /no DATA9 record/);
});
