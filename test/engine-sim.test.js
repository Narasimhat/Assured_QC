// End-to-end: synthetic control/edited pairs with known allele mixtures, analysed against the design that
// produced them. Tolerances are stated points of percentage difference from the known mixture.
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { parseAbif } from "../src/abif.js";
import { analysePair } from "../src/analyze.js";
import { loadSpec, simDir, loadTruth } from "./helpers.js";

const truth = loadTruth();
const oracle = JSON.parse(readFileSync(path.join(simDir, "oracle_ice.json"), "utf8")).results;
const traces = new Map();
const trace = (name) => { if (!traces.has(name)) traces.set(name, parseAbif(readFileSync(path.join(simDir, name)), name)); return traces.get(name); };
const specFor = (name) => loadSpec(name.startsWith("nr2f2") ? "nr2f2_2xha_ssodn" : "apoe_r176c_snp");
const results = new Map();
const run = (name) => {
  if (!results.has(name)) {
    const entry = truth[name];
    results.set(name, analysePair({ control: trace(entry.control), edited: trace(entry.edited), spec: specFor(name) }));
  }
  const result = results.get(name);
  assert.ok(result.ok, result.error);
  return result;
};
const frac = (name, label) => 100 * (truth[name].alleles.find((a) => a.label === label)?.fraction || 0);
const near = (actual, expected, tolerance, message) => assert.ok(Math.abs(actual - expected) <= tolerance, `${message}: ${actual} vs ${expected} (+/-${tolerance})`);

test("knockout pool: editing, wild type and knockout score match the known mixture and ICE", () => {
  const { summary } = run("nr2f2_ko_pool_g1");
  near(summary.wtPct, 25, 4, "wild type");
  near(summary.editedPct, 75, 4, "edited");
  near(summary.koScorePct, 53, 4, "KO score (frameshift or >=21 bp)");
  near(summary.editedPct, oracle.nr2f2_ko_pool_g1.summary.ice, 4, "vs ICE");
  assert.equal(summary.intendedEditPct, 0);
});

test("knockout clones: heterozygous, compound heterozygous and homozygous", () => {
  near(run("nr2f2_ko_clone_het_wt_m1").summary.wtPct, 50, 5, "het wild type");
  near(run("nr2f2_ko_clone_compound").summary.editedPct, 100, 3, "compound het edited");
  near(run("nr2f2_ko_clone_hom_m1").summary.editedPct, 100, 3, "hom edited");
  const compound = run("nr2f2_ko_clone_compound").decomposition.contributions.slice(0, 2).map((c) => [c.size, Math.round(c.fraction * 10) / 10]);
  assert.deepEqual(compound.map(([size]) => size).sort((a, b) => a - b), [-7, 1]);
  const hom = run("nr2f2_ko_clone_hom_m1").decomposition.contributions[0];
  assert.equal(hom.size, -1);
  assert.ok(hom.fraction > 0.93);
});

test("two guides: the deletion between the cuts is found, with its size", () => {
  const result = run("nr2f2_del_two_guides");
  near(result.summary.betweenCutsPct, 25, 5, "deletion between cuts");
  const top = result.decomposition.contributions.find((c) => c.kind === "deletion_between_cuts");
  assert.ok(top && Math.abs(top.size) >= 26 && Math.abs(top.size) <= 32, `size ${top?.size}`);
  near(result.summary.editedPct, 70, 5, "edited");
});

test("small-tag ssODN (96 bp insert): intended edit in a pool and in clones", () => {
  near(run("nr2f2_ha_hdr_pool").summary.intendedEditPct, 30, 4, "pool intended edit");
  near(run("nr2f2_ha_hdr_pool").summary.intendedEditPct, oracle.nr2f2_ha_hdr_pool.summary.hdr_pct, 4, "vs ICE HDR");
  near(run("nr2f2_ha_hdr_clone_het").summary.intendedEditPct, 50, 5, "het clone");
  near(run("nr2f2_ha_hdr_clone_hom").summary.intendedEditPct, 100, 6, "hom clone");
});

test("SNP with silent blocking changes: partial conversion is separated from the intended edit", () => {
  const pool = run("apoe_snp_hdr_pool").summary;
  near(pool.intendedEditPct, 35, 6, "pool intended edit");
  near(pool.partialConversionPct, 10, 6, "pool blocking-only");
  near(pool.wtPct, 20, 6, "pool wild type");
  const het = run("apoe_snp_clone_het").summary;
  near(het.intendedEditPct, 50, 6, "het clone intended");
  near(het.wtPct, 50, 6, "het clone wild type");
  const partial = run("apoe_snp_clone_partial_het").summary;
  near(partial.intendedEditPct, 50, 7, "clone with blocking changes only on one allele: intended edit");
  near(partial.partialConversionPct, 50, 7, "clone with blocking changes only on one allele: partial");
  // ICE counts the blocking-only allele as knock-in
  assert.ok(oracle.apoe_snp_clone_partial_het.summary.hdr_pct > 90);
  assert.ok(partial.intendedEditPct < 60);
});

test("the fit explains the data and the abundances sum to about one", () => {
  for (const name of Object.keys(truth)) {
    const { decomposition } = run(name);
    if (decomposition.ssWt > 0.5) assert.ok(decomposition.r2 > (truth[name].kind.startsWith("noisy") ? 0.8 : 0.9), `${name} R2 ${decomposition.r2}`); // R2 is relative to the wild-type-only fit, so it means little when nothing is edited
    assert.ok(Math.abs(decomposition.rawSum - 1) < 0.1, `${name} sum ${decomposition.rawSum}`);
  }
});

test("a pair of identical traces shows no editing", () => {
  const control = trace("nr2f2_control.ab1");
  const result = analysePair({ control, edited: control, spec: specFor("nr2f2") });
  assert.ok(result.ok);
  near(result.summary.wtPct, 100, 2, "wild type");
  assert.ok(result.summary.editedPct < 2.5, `spurious editing ${result.summary.editedPct}%`);
});
