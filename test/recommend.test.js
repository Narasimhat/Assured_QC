import test from "node:test";
import assert from "node:assert/strict";
import { mkdtempSync, copyFileSync, readFileSync, existsSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { main as cli } from "../cli/assured-qc.mjs";
import { alleleRisk, supportOf, recommendClones, withRecommendations, panelOf, recommendationCsv, applyObservations, RISK_DEFAULTS } from "../src/recommend.js";
import { buildPlateReportHtml } from "../src/plateReport.js";
import { buildResultsCsv } from "../src/report.js";

const here = path.dirname(fileURLToPath(import.meta.url)); const root = path.join(here, "..");

// allele rows as the functional module writes them
const NMD = (label, fraction = 1, extra = {}) => ({ label, fraction, class: "lof_nmd", lof: true, details: { nmd: { predicted: true, borderline: false }, nearSplice: false, ...extra } });
const NMD_BORDER = (label, fraction = 1) => ({ label, fraction, class: "lof_nmd", lof: true, details: { nmd: { predicted: true, borderline: true }, nearSplice: true } });
const SPLICE = (label, fraction = 1) => ({ label, fraction, class: "lof_splice", lof: true, details: { splice: { skipping: { class: "lof_truncation", lof: true }, retention: { class: "lof_nmd", lof: true } } } });
const INFRAME = (label, fraction = 1) => ({ label, fraction, class: "inframe", lof: false, details: {} });
const WT = (fraction) => ({ label: "Wild type", fraction, class: "wt", lof: false });

// a knockout clone with its functional rows: alleles = [[row, ...]] shares come from the rows' fractions
const clone = (name, rows, over = {}) => {
  const alleles = rows.filter((r) => r.class !== "wt").map((r) => ({ label: r.label, kind: "deletion_between_cuts", size: -1 }));
  return { ok: true, kind: "trace-pair", name, workflow: "deletion", sampleType: "clone", confidence: { tier: "high" }, warnings: [], reads: undefined,
    summary: { wtPct: rows.some((r) => r.class === "wt") ? 50 : 0, editedPct: 100, koScorePct: 100, functionalKoPct: 100 }, intervals: { editedPct: [96, 100], wtPct: [0, 3] },
    genotype: { category: alleles.length > 1 ? "compound_heterozygous" : "homozygous_indel", summary: alleles.map((a) => a.label).join(" + "), alleles },
    functional: { call: { status: "predicted_knockout_with_caveats", text: "x" }, caveats: ["c"], confirm: [], rows },
    decision: { decision: "hold", reasons: ["Homozygous", "caveats"], next: "check protein" }, ...over };
};

test("allele risk follows the mechanism: clean NMD < borderline NMD < truncation < splice-dependent < undecided < wild type", () => {
  assert.equal(alleleRisk(NMD("a")).risk, RISK_DEFAULTS.lof_nmd);
  assert.equal(alleleRisk(NMD_BORDER("a")).risk, RISK_DEFAULTS.lof_nmd_borderline);
  assert.ok(Math.abs(alleleRisk(SPLICE("a")).risk - (RISK_DEFAULTS.lof_truncation + RISK_DEFAULTS.crypticSplice)) < 1e-9, "the weaker of the two splice outcomes plus the cryptic-site margin");
  assert.equal(alleleRisk(INFRAME("a")).risk, RISK_DEFAULTS.not_lof);
  assert.equal(alleleRisk(WT(1)).risk, 1);
  assert.equal(alleleRisk({ class: "truncation_uncertain", lof: false, label: "x" }).risk, RISK_DEFAULTS.uncertain);
  const order = [NMD("a"), NMD_BORDER("a"), { label: "t", class: "lof_truncation", lof: true }, SPLICE("a"), { label: "u", class: "frameshift_no_stop", lof: false }, INFRAME("a"), WT(1)].map((r) => alleleRisk(r).risk);
  assert.deepEqual([...order].sort((a, b) => a - b), order);
  assert.equal(alleleRisk(SPLICE("a"), { crypticSplice: 0 }).risk, RISK_DEFAULTS.lof_truncation, "weights can be replaced once the calibration is known");
});

test("support is 100 x (1 - share-weighted risk): nothing short of a protein result reaches 100", () => {
  assert.equal(supportOf(clone("c", [NMD("a")])).score, 85);
  assert.equal(supportOf(clone("c", [SPLICE("a")])).score, 50);
  assert.equal(supportOf(clone("c", [NMD_BORDER("a", 0.5), SPLICE("b", 0.5)])).score, 60);
  assert.equal(supportOf(clone("c", [NMD("a", 0.5), WT(0.5)])).tier, "unsupported");
  assert.equal(supportOf(clone("c", [NMD("a")])).tier, "firm");
  assert.equal(supportOf(clone("c", [SPLICE("a")])).tier, "possible");
  assert.equal(supportOf(clone("c", [NMD_BORDER("a")])).tier, "likely");
  assert.ok(supportOf(clone("c", [NMD("a", 0.9), NMD("tiny", 0.02)])).score === 85, "an allele under 5% is ignored");
});

test("ranking: support first, then homozygous before compound, then trace quality; one clone per genotype goes to the panel", () => {
  const results = [
    clone("c1", [SPLICE("del79")]), clone("c2", [SPLICE("del79")]), clone("c3", [NMD_BORDER("del77", 0.5), SPLICE("del79", 0.5)], { warnings: ["The edited trace has low quality scores over much of its length; treat the result with caution."] }),
    clone("c4", [NMD("del17")]), clone("c5", [WT(0.9), NMD("x", 0.1)], { decision: { decision: "reject", reasons: ["wt"], next: "discard" } }),
    clone("c6", [NMD("del17")], { decision: { decision: "re-sequence", reasons: ["noisy"], next: "again" } }),
    clone("c7", [SPLICE("del79", 0.5), SPLICE("del84", 0.5)]),
  ];
  const rec = recommendClones(results); const by = Object.fromEntries(results.map((r, i) => [r.name, rec[i]]));
  assert.equal(by.c4.rank, 1); assert.equal(by.c3.rank, 2); assert.equal(by.c3.score, 60);
  assert.equal(by.c1.rank, 3, "homozygous c1 outranks compound c7 at the same support"); assert.equal(by.c2.rank, 4); assert.equal(by.c7.rank, 5);
  assert.equal(by.c4.priority, "first"); assert.equal(by.c3.priority, "first"); assert.equal(by.c1.priority, "first");
  assert.equal(by.c2.priority, "duplicate"); assert.equal(by.c2.duplicateOf, "c1"); assert.match(by.c2.action, /backup/i); assert.equal(by.c1.groupOf, 2);
  assert.equal(by.c7.priority, "reserve", "the panel holds three genotypes");
  assert.equal(by.c5.priority, "discard"); assert.equal(by.c5.rank, null); assert.equal(by.c6.priority, "resequence"); assert.equal(by.c6.rank, null);
  assert.ok(by.c3.flags.some((f) => /trace quality/.test(f)) && by.c3.flags.some((f) => /one read/.test(f)));
  const dup = recommendClones([clone("d", [SPLICE("del79", 0.25), SPLICE("del79", 0.25), NMD("del17", 0.5)])])[0]; assert.equal(dup.why.length, 2); assert.match(dup.why[0] + dup.why[1], /50% of the signal/);
  assert.deepEqual(panelOf(withRecommendations(results)), ["c4", "c3", "c1"]);
  assert.deepEqual(panelOf(withRecommendations(results, { panel: 1 })), ["c4"]);
});

test("equal support: the clone without a trace-quality warning goes first, and ties fall back to the sample name in natural order", () => {
  const low = ["The edited trace has low quality scores over much of its length; treat the result with caution."];
  const results = [clone("clone10", [SPLICE("a")]), clone("clone2", [SPLICE("b")], { warnings: low }), clone("clone3", [SPLICE("c")])];
  const rec = recommendClones(results);
  assert.deepEqual(results.map((r, i) => [r.name, rec[i].rank]).sort((a, b) => a[1] - b[1]).map((x) => x[0]), ["clone3", "clone10", "clone2"]);
});

test("without a gene record the support is the frameshift proxy and says so; edit workflows rank by the intended-edit lower bound", () => {
  const ko = clone("k", [NMD("a")]); delete ko.functional; ko.summary = { wtPct: 5, editedPct: 95, koScorePct: 100 };
  const r = recommendClones([ko])[0]; assert.equal(r.score, 95); assert.match(r.basis, /^KO score/); assert.match(r.basis, /no gene record/);
  const snp = { ok: true, kind: "trace-pair", name: "s", workflow: "snp", sampleType: "clone", confidence: { tier: "high" }, warnings: [], summary: { intendedEditPct: 90 }, intervals: { intendedEditPct: [86, 94] }, genotype: { category: "homozygous_edit", alleles: [{ label: "Intended edit", kind: "edit", size: 0 }] }, decision: { decision: "accept", reasons: [], next: "" } };
  const s = recommendClones([snp])[0]; assert.equal(s.score, 86); assert.equal(s.tier, "firm"); assert.deepEqual(s.testWith, []);
  assert.deepEqual(recommendClones([{ ok: true, kind: "junction", name: "j", decision: { decision: "accept", reasons: [] } }, { ok: false, name: "e", decision: { decision: "re-sequence", reasons: ["x"] } }]).map((x) => x.priority), ["none", "none"]);
});

test("protein results: absent confirms a clean knockout, present rejects, truncated and reduced hold; the earlier decision is kept", () => {
  const results = [clone("a", [SPLICE("d79")]), clone("b", [SPLICE("d78")]), clone("c", [SPLICE("d77")]), clone("d", [SPLICE("d84")]), clone("e", [NMD("d17")], { workflow: "snp" }), clone("f", [SPLICE("d1")]),
    clone("g", [SPLICE("d2")], { summary: { wtPct: 30, editedPct: 70, koScorePct: 100 }, intervals: { editedPct: [60, 80], wtPct: [20, 40] } })];
  const applied = applyObservations(results, [{ sample: "a", protein: "absent", notes: "N-terminal antibody" }, { sample: "b", protein: "present" }, { sample: "c", protein: "truncated" }, { sample: "d", protein: "reduced" }, { sample: "e", protein: "absent" }, { sample: "nothere", protein: "absent" }, { sample: "f", protein: "maybe" }, { sample: "g", protein: "absent" }]);
  const by = Object.fromEntries(applied.results.map((r) => [r.name, r]));
  assert.equal(by.a.decision.decision, "accept"); assert.equal(by.a.decision.before.decision, "hold"); assert.match(by.a.decision.reasons.join(" "), /Protein absent.*N-terminal antibody/);
  assert.equal(by.b.decision.decision, "reject"); assert.equal(by.c.decision.decision, "hold"); assert.match(by.c.decision.reasons.join(" "), /shorter product/); assert.equal(by.d.decision.decision, "hold");
  assert.equal(by.e.decision.decision, "hold", "protein results change knockout and deletion samples only"); assert.equal(by.e.confirmation.protein, "absent");
  assert.equal(by.g.decision.decision, "hold", "absent protein does not accept a clone that still carries wild type in the sequence"); assert.match(by.g.decision.reasons.join(" "), /does not show a clean knockout/);
  assert.equal(by.f.decision.decision, "hold"); assert.equal(applied.matched, 6);
  assert.ok(applied.problems.some((p) => /nothere/.test(p)) && applied.problems.some((p) => /"maybe"/.test(p)));
  assert.equal(results[0].decision.decision, "hold", "the inputs are not modified");
  const rec = recommendClones(applied.results); const a = rec[0];
  assert.equal(a.score, 100); assert.equal(a.tier, "confirmed"); assert.equal(a.priority, "confirmed"); assert.equal(a.rank, 1);
  assert.equal(rec[1].priority, "discard", "a clone whose protein is present is not carried forward");
  const rec2 = recommendClones(applied.results.filter((r) => ["a", "c", "d"].includes(r.name))); const byc = Object.fromEntries(applied.results.filter((r) => ["a", "c", "d"].includes(r.name)).map((r, i) => [r.name, rec2[i]]));
  assert.equal(byc.a.priority, "confirmed"); assert.equal(byc.c.priority, "caution", "a measured truncated product is not carried forward"); assert.equal(byc.d.priority, "caution"); assert.ok(byc.c.score <= 20 && byc.c.rank > byc.a.rank);
});

test("recommendation CSV: ranked, quoted, one row per clone with a rank or a verdict", () => {
  const results = withRecommendations([clone("b", [SPLICE("d, 79")]), clone("a", [NMD("d17")]), clone("x", [NMD("d17")], { decision: { decision: "reject", reasons: ["wt"], next: "" } })]);
  const lines = recommendationCsv(results).trim().split("\n"); assert.equal(lines.length, 4);
  assert.match(lines[0], /^rank,sample,well,decision,priority,knockout_support,support_tier/); assert.match(lines[1], /^1,a,/); assert.match(lines[2], /^2,b,/); assert.match(lines[2], /"[^"]*d, 79[^"]*"/); assert.match(lines[3], /^,x,.*,discard,/);
  const csv = buildResultsCsv("deletion", results).split("\n"); assert.ok(/knockout_support,priority_rank,recommendation$/.test(csv[0])); assert.ok(csv[1].split(",").length >= csv[0].split(",").length);
  const html = buildPlateReportHtml({ results }); assert.match(html, /Order to carry forward and confirm \(2\)/); assert.match(html, /First choice: a, b/); assert.ok(!/<script/i.test(html));
  assert.doesNotMatch(buildPlateReportHtml({ results: results.map(({ recommendation, ...r }) => r) }), /Order to carry forward/);
});

test("command line: run ranks the clones and writes recommendation.csv; recommend re-ranks a finished run with protein results", async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "qc-rec-")); const sim = (n) => path.join(root, "fixtures", "sim", n);
  const files = { "apoe_WT.700.ab1": "apoe_control.ab1", "apoe_A01_hom.700.ab1": "apoe_snp_clone_hom_edited.ab1", "apoe_A02_het.700.ab1": "apoe_snp_clone_het_edited.ab1" };
  Object.entries(files).forEach(([to, from]) => copyFileSync(sim(from), path.join(dir, to)));
  const out = path.join(dir, "out"); const logged = []; const original = console.log; console.log = (...a) => logged.push(a.join(" "));
  try {
    assert.equal(await cli(["run", dir, "--design", path.join(root, "fixtures", "design", "apoe_r176c_snp.design.json"), "--threads", "2", "--out", out]), 0);
    assert.ok(existsSync(path.join(out, "recommendation.csv")));
    const summary = readFileSync(path.join(out, "summary.svg"), "utf8"); assert.match(summary, /^<svg /); assert.match(summary, /apoe_A01_hom\.700/); assert.match(summary, /#1 carry forward/);
    assert.match(readFileSync(path.join(out, "report.html"), "utf8"), /<h2>Summary<\/h2><div style="overflow-x:auto"><svg /);
    const results = JSON.parse(readFileSync(path.join(out, "results.json"), "utf8")); assert.ok(results.every((r) => r.recommendation && r.recommendation.priority !== undefined));
    const by = Object.fromEntries(results.map((r) => [r.name, r])); assert.equal(by["apoe_A01_hom.700"].recommendation.rank, 1); assert.ok(by["apoe_A02_het.700"].recommendation.rank > 1);
    assert.match(logged.join("\n"), /Take forward first: apoe_A01_hom\.700/);
    // recommend on the saved run, with a protein table (snp samples are not changed by protein results)
    writeFileSync(path.join(dir, "protein.csv"), "sample,protein\napoe_A02_het.700,absent\n"); const target = path.join(dir, "rec.csv");
    assert.equal(await cli(["recommend", path.join(out, "results.json"), "--observations", path.join(dir, "protein.csv"), "--panel", "1", "--out", target]), 0);
    assert.ok(existsSync(target)); const rec = JSON.parse(readFileSync(target.replace(/\.csv$/, ".json"), "utf8")); assert.equal(rec.length, 2);
    assert.equal(await cli(["recommend", path.join(dir, "missing.json")]), 1);
  } finally { console.log = original; }
});

test("a wanted edit: clones of the same genotype are replicates and all count; clones that do not carry the edit are not recommended and do not take a rank", () => {
  const snp = (name, lower, decision = "accept", over = {}) => ({ ok: true, kind: "trace-pair", name, workflow: "snp", sampleType: "clone", confidence: { tier: "high" }, warnings: [], summary: { intendedEditPct: lower + 2, editedPct: 100 }, intervals: { intendedEditPct: [lower, 100] },
    genotype: { category: "homozygous_edit", alleles: [{ label: "Intended edit", kind: "edit", size: 0 }] }, decision: { decision, reasons: [], next: "" }, ...over });
  const results = [snp("a", 98), snp("b", 97), snp("c", 96), snp("d", 84, "hold"), snp("e", 0, "hold", { genotype: { category: "homozygous_indel", alleles: [{ label: "-13 bp", kind: "indel", size: -13 }] } }), snp("f", 10, "hold")];
  const rec = recommendClones(results, { panel: 3 }); const by = Object.fromEntries(results.map((r, i) => [r.name, rec[i]]));
  assert.deepEqual(["a", "b", "c"].map((n) => by[n].priority), ["first", "first", "first"], "three clones of one genotype are all first choices");
  assert.equal(by.d.priority, "reserve"); assert.match(by.d.action, /wanted edit/); assert.match(by.a.action, /wanted edit/);
  assert.equal(by.e.priority, "not_recommended"); assert.equal(by.e.rank, null); assert.equal(by.f.priority, "not_recommended"); assert.match(by.e.action, /does not support the wanted edit/);
  assert.deepEqual(["a", "b", "c", "d"].map((n) => by[n].rank), [1, 2, 3, 4], "ranks are contiguous: not recommended clones take none");
  assert.deepEqual(panelOf(withRecommendations(results, { panel: 3 })), ["a", "b", "c"]);
  const ko = clone("k", [NMD("x", 0.5), WT(0.5)]); assert.equal(recommendClones([ko])[0].priority, "not_recommended", "a knockout clone with wild type is not recommended");
});

test("a knockout plate without a gene record: the default panel is three (not zero), the KO score is named, up to two clones of one deletion count, ties go to the narrower interval", () => {
  const ko = (name, alleles, lo, over = {}) => ({ ok: true, kind: "trace-pair", name, workflow: "deletion", sampleType: "clone", confidence: { tier: "high" }, warnings: [], summary: { wtPct: 0, editedPct: 100, koScorePct: 100 }, intervals: { editedPct: [98, 100], koScorePct: [lo, 100] },
    genotype: { category: alleles.length === 1 ? "homozygous_indel" : "compound_heterozygous", alleles: alleles.map((a) => ({ label: `Deletion between cuts (${a} bp)`, kind: "deletion_between_cuts", size: -a })) },
    contributions: alleles.map((a) => ({ label: `Deletion between cuts (${a} bp)`, kind: "deletion_between_cuts", size: -a, fraction: 1 / alleles.length })), decision: { decision: "accept", reasons: [], next: "" }, ...over });
  const results = [ko("a1", [79], 98.8), ko("a2", [78], 99), ko("a3", [79], 99), ko("a4", [79], 99), ko("a5", [79], 99), ko("c1", [84, 79], 98.3), ko("c2", [77, 79], 86.5)];
  const rec = withRecommendations(results); const by = Object.fromEntries(rec.map((r) => [r.name, r.recommendation]));
  assert.deepEqual(panelOf(rec), ["a2", "a3", "a4"], "no panel option given: three clones, the narrower intervals first, at most two of the same deletion");
  assert.match(by.a2.basis, /^KO score/); assert.equal(by.a5.priority, "duplicate"); assert.equal(by.a1.priority, "duplicate"); assert.equal(by.c1.priority, "reserve");
  assert.match(by.a4.action, /second clone with the same deletion/); assert.deepEqual(panelOf(withRecommendations(results, { panel: undefined })), ["a2", "a3", "a4"], "an undefined panel option is not a zero");
  assert.deepEqual(panelOf(withRecommendations(results, { panel: 1 })), ["a2"]);
});
