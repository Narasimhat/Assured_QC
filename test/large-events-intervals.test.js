import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { analysePair } from "../src/analyze.js";
import { specFromControl } from "../src/manualSpec.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const here = path.dirname(fileURLToPath(import.meta.url));
const extended = JSON.parse(readFileSync(path.join(here, "..", "fixtures", "sim", "truth_extended.json"), "utf8"));
const base = loadTruth();
const specFor = (name) => loadSpec(/nr2f2/.test(name) ? "nr2f2_2xha_ssodn" : "apoe_r176c_snp");
const run = (entry, name, extra = {}) => analysePair({ control: loadTrace(entry.control), edited: loadTrace(entry.edited), spec: specFor(name), options: { sampleType: "clone", ...extra } });
const share = (result, size) => result.decomposition.contributions.filter((c) => c.size === size && (c.kind === "indel")).reduce((s, c) => s + c.fraction, 0);

test("a 60-base deletion at one cut is found by the adaptive step, not listed in the base library", () => {
  const entry = extended.ext_nr2f2_large_del60_het; const result = run(entry, "nr2f2");
  assert.ok(result.ok, result.error);
  assert.ok(share(result, -60) > 0.4 && share(result, -60) < 0.6, `-60 share ${share(result, -60)}`);
  assert.ok(Math.abs(result.summary.editedPct - 50) < 6, `edited ${result.summary.editedPct}`);
  assert.ok(result.decomposition.pursuit, "the adaptive step added the allele");
  const off = run(entry, "nr2f2", { pursuit: false });
  assert.ok(share(off, -60) === 0, "without the adaptive step the allele is not in the library");
});

test("a 6-base insertion of unknown sequence is recognised by its shift", () => {
  const result = run(extended.ext_nr2f2_ins6_het, "nr2f2"); assert.ok(result.ok, result.error);
  assert.ok(share(result, 6) > 0.4 && share(result, 6) < 0.6, `+6 share ${share(result, 6)}`);
  assert.ok(Math.abs(result.summary.editedPct - 50) < 6);
});

test("wild type: no false large allele, and the interval says how much could be hidden", () => {
  const result = run(extended.ext_nr2f2_wt_replicate, "nr2f2"); assert.ok(result.ok, result.error);
  assert.ok(result.summary.editedPct < 2, `edited ${result.summary.editedPct}`);
  assert.equal(result.decomposition.pursuit, undefined);
  assert.equal(result.intervals.detected, false);
  // the profile is evaluated on a 2.5-point grid and the interval has a 1-point floor, so a clean wild type reports at most 6 points that could be hidden
  assert.ok(result.intervals.hiddenUpToPct <= 6, `could hide up to ${result.intervals.hiddenUpToPct}`);
});

test("a marker upstream of an indel on the same allele cannot be phased: the interval spans both readings", () => {
  const entry = extended.ext_apoe_phase_ambiguous; const control = loadTrace(entry.control);
  const built = specFromControl({ control, guides: entry.guides.join(", "), donor: entry.donor, gene: "AMB" }); assert.ok(!built.error, built.error);
  const result = analysePair({ control, edited: loadTrace(entry.edited), spec: built.spec, options: { sampleType: "clone" } }); assert.ok(result.ok, result.error);
  // truth: 50% wild type + 50% (SNP and +1). The same bases read as 50% SNP-only + 50% +1, which has 100% edited.
  const [lo, hi] = result.intervals.editedPct;
  assert.ok(lo <= 52 && hi >= 90, `edited interval ${lo}-${hi} must cover the truth (50) and the alternative reading (near 100)`);
  const [ilo, ihi] = result.intervals.intendedEditPct;
  assert.ok(ilo <= 2 && ihi >= 40, `intended-edit interval ${ilo}-${ihi}`);
});

test("when markers flank the indel the point estimate reads one edited allele, and the interval still admits the partial-conversion reading", () => {
  const entry = extended.ext_apoe_flanked; const result = run(entry, "apoe"); assert.ok(result.ok, result.error);
  const [lo, hi] = result.intervals.editedPct;
  // 50% wild type + 50% (markers and indel) gives the same base composition as 50% (left marker only) + 50% (indel and right marker), so flanking markers
  // do not make the phase identifiable from one trace; the interval must contain the truth and may extend to the partial-conversion reading
  assert.ok(lo <= 50 && hi >= 50, `edited interval ${lo}-${hi}`);
  assert.ok(Math.abs(result.summary.editedPct - 50) < 6, `edited ${result.summary.editedPct}`);
});

test("intervals contain the truth for the existing clone fixtures", () => {
  for (const name of ["nr2f2_ko_clone_het_wt_m1", "nr2f2_ko_clone_hom_m1", "nr2f2_ko_clone_compound", "nr2f2_ha_hdr_clone_het", "nr2f2_ha_hdr_clone_hom", "apoe_snp_clone_het", "apoe_snp_clone_hom"]) {
    const entry = base[name]; const result = run(entry, name); assert.ok(result.ok, result.error);
    const wt = entry.alleles.filter((a) => a.label === "wt").reduce((s, a) => s + a.fraction, 0) * 100; const edited = 100 - wt;
    const [lo, hi] = result.intervals.editedPct;
    assert.ok(lo - 0.1 <= edited && edited <= hi + 0.1, `${name}: truth ${edited} not in ${lo}-${hi}`);
    // a tag design has partial-conversion alleles that the composition cannot rule out (50% wild type + 50% edit reads like 50% insert-only + 50% markers-only), so its interval may be wide
    if (!/ha_hdr/.test(name)) assert.ok(hi - lo < 30, `${name}: interval ${lo}-${hi} is too wide for a clean trace`);
  }
});

test("a donor edit that also carries a small indel is a library allele (flanked case is read as one allele)", () => {
  const result = run(extended.ext_apoe_flanked, "apoe");
  assert.ok(result.decomposition.contributions.some((c) => c.kind === "edit_indel" && c.fraction > 0.35), JSON.stringify(result.decomposition.contributions.slice(0, 3).map((c) => [c.label, c.kind, c.fraction])));
});
