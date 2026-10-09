// Regressions found when the engine was run on 172 control/edited pairs from laboratory traces (docs/validation.md).
import test from "node:test";
import assert from "node:assert/strict";
import { analysePair, summarise } from "../src/analyze.js";
import { loadSpec, loadTrace } from "./helpers.js";

const spec = loadSpec("nr2f2_2xha_ssodn");
const control = loadTrace("nr2f2_control.ab1");
const het = loadTrace("nr2f2_ha_hdr_clone_het_edited.ab1");

// Base calls replaced by pseudo-random bases over the first `count` positions, as in a read whose first few hundred calls
// are poor: the first stretch of the control's good region then cannot be aligned to the edited trace.
const scrambledStart = (trace, count) => {
  let seed = 11; const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const calls = trace.calls.split("");
  for (let i = 0; i < Math.min(count, calls.length); i += 1) calls[i] = "ACGT"[Math.floor(rnd() * 4)];
  return { ...trace, calls: calls.join("") };
};

test("an edited read whose start cannot be called is still aligned from the part that matches", () => {
  const edited = scrambledStart(het, 150);
  const result = analysePair({ control, edited, spec, options: { sampleType: "clone" } });
  assert.ok(result.ok, result.error);
  assert.ok(result.warnings.some((w) => /bases of the control's upstream stretch align/.test(w)), result.warnings.join(" | "));
  assert.ok(Math.abs(result.summary.intendedEditPct - 50) < 8, `intended ${result.summary.intendedEditPct}`);
});

test("an unrelated edited trace is still refused", () => {
  const edited = scrambledStart(het, het.calls.length);
  const result = analysePair({ control, edited, spec, options: { sampleType: "clone" } });
  assert.equal(result.ok, false);
  assert.match(result.error, /does not align to the control upstream of the cut/);
});

const decomposition = (noiseFraction) => ({
  noiseFraction,
  contributions: [
    { kind: "wt", size: 0, fraction: 0.96 },
    { kind: "indel", size: -2, fraction: 0.03 },
    { kind: "indel", size: 3, fraction: 0.01 },
  ],
});

test("a clean fit is summarised as before", () => {
  const s = summarise(decomposition(0.04));
  assert.equal(s.wtPct, 96);
  assert.equal(s.editedPct, 4);
  assert.equal(s.unexplainedPct, 0);
});

test("signal that no allele explains is counted as not wild type, not renormalised away", () => {
  // A complex pool: the wild type explains 6% of the signal, hundreds of rare alleles are absorbed as background.
  const s = summarise(decomposition(0.94));
  assert.equal(s.unexplainedPct, 84);
  assert.ok(s.wtPct < 16, `wild type ${s.wtPct}`);
  assert.ok(s.editedPct > 84, `edited ${s.editedPct}`);
  assert.ok(Math.abs(s.wtPct + s.indelPct + s.unexplainedPct - 100) < 0.5);
});
