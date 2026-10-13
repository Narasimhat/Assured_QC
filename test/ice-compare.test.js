import test from "node:test";
import assert from "node:assert/strict";
import { mkdtempSync, mkdirSync, writeFileSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { appScore, appAlleles, iceJobsFor, compareWithIce, compareRun, iceSummary, iceCsv, iceScatterSvg, iceSectionHtml, iceRawRecords, matchIceRaw, parseIceJson, ICE_FIT_OK, ICE_DIFFERENCE } from "../src/iceCompare.js";
import { checkIce, runIceJobs, ICE_SCRIPT } from "../cli/ice-runner.mjs";
import { convertSvg } from "../src/svgPdf.js";
import { buildQcReportHtml } from "../src/report.js";
import { analyseSample } from "../src/batch.js";
import { withRecommendations } from "../src/recommend.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

// simulated results: a knockout clone with two deletions that differ by one base, and a point edit
const ko = (name, over = {}) => ({ ok: true, name, workflow: "deletion", sampleType: "clone", summary: { koScorePct: 100, wtPct: 0 }, intervals: { koScorePct: [97, 100] },
  contributions: [{ kind: "deletion_between_cuts", size: -77, fraction: 0.5 }, { kind: "deletion_between_cuts", size: -78, fraction: 0.5 }], ...over });
const edit = (name, pct, over = {}) => ({ ok: true, name, workflow: "snp", sampleType: "clone", summary: { intendedEditPct: pct, wtPct: 100 - pct }, intervals: { intendedEditPct: [pct - 5, 100] }, contributions: [], ...over });
const run = (ko_score, rsq, contribs, extra = {}) => ({ ok: true, donorName: "", summary: { ice: ko_score, ko_score, rsq }, contribs: contribs.map(([indel, rel]) => ({ indel, rel })), ...extra });
const donorRun = (donorName, hdr_pct, rsq = 0.98) => ({ ok: true, donorName, summary: { ice: hdr_pct, hdr_pct, rsq }, contribs: [] });

test("the app score is the KO score for knockouts and the intended-edit share for edits; same-size alleles are merged", () => {
  assert.deepEqual(appScore(ko("a")), { value: 100, interval: [97, 100], label: "KO score" });
  assert.equal(appScore(ko("a", { summary: { koScorePct: 100, wtPct: 25 } })).value, 75);
  assert.deepEqual(appScore(edit("e", 88)), { value: 88, interval: [83, 100], label: "HDR score" });
  const slid = ko("s", { contributions: [{ kind: "deletion_between_cuts", size: -79, fraction: 0.41 }, { kind: "deletion_between_cuts", size: -79, fraction: 0.08 }, { kind: "deletion_between_cuts", size: -84, fraction: 0.49 }, { kind: "wt", size: 0, fraction: 0.02 }] });
  assert.equal(appAlleles(slid).length, 2, "one entry per size"); assert.ok(appAlleles(slid).some((a) => a.size === -84)); assert.ok(Math.abs(appAlleles(slid).find((a) => a.size === -79).share - 0.49) < 1e-9);
});

test("agreement: scores within 15 points and a sound fit", () => {
  const row = compareWithIce(ko("a", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] }), [run(99, 0.99, [[-78, 0.99]])]);
  assert.equal(row.status, "agree"); assert.equal(row.difference, 1); assert.deepEqual(row.missing, []); assert.equal(row.note, "");
});

test("poor ICE fit: not a disagreement, and the allele ICE does not place is named (compound heterozygote)", () => {
  const row = compareWithIce(ko("het"), [run(58, 0.63, [[-78, 0.41], [5, 0.05], [4, 0.03]])]);
  assert.equal(row.status, "ice_fit_poor"); assert.equal(row.app, 100); assert.equal(row.ice, 58); assert.equal(row.iceR2, 0.63);
  assert.deepEqual(row.missing.map((a) => a.size), [-77]);
  assert.match(row.note, /not comparable/); assert.match(row.note, /-77 \(50%\), which ICE does not/); assert.ok(ICE_FIT_OK === 0.8);
});

test("a real difference with a sound fit is flagged for a look at the trace; a small one is not", () => {
  const big = compareWithIce(ko("x", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] }), [run(60, 0.95, [[-78, 0.6]])]);
  assert.equal(big.status, "differ"); assert.equal(big.difference, 40); assert.match(big.note, /look at the trace/);
  const edge = compareWithIce(ko("y", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] }), [run(100 - ICE_DIFFERENCE + 1, 0.95, [[-78, 0.9]])]);
  assert.equal(edge.status, "agree");
  assert.equal(compareWithIce(ko("z", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] }), [run(100 - ICE_DIFFERENCE, 0.95, [[-78, 0.85]])]).status, "differ");
});

test("point edits: one ICE run per donor, the best donor is compared, and the others are listed", () => {
  const row = compareWithIce(edit("e", 90), [donorRun("ssODN1", 45), donorRun("ssODN2", 88)]);
  assert.equal(row.ice, 88); assert.equal(row.donor, "ssODN2"); assert.equal(row.status, "agree"); assert.match(row.note, /Best of 2 donors.*ssODN1: 45, ssODN2: 88/);
});

test("failures and samples ICE was not run for stay visible", () => {
  assert.equal(compareWithIce(ko("f"), [{ ok: false, donorName: "", error: "ValueError: no cut site" }]).status, "ice_failed");
  assert.match(compareWithIce(ko("f"), [{ ok: false, error: "ValueError: no cut site" }]).note, /no cut site/);
  const none = compareWithIce(ko("n"), [], "ICE supports SpCas9 only"); assert.equal(none.status, "not_run"); assert.equal(none.note, "ICE supports SpCas9 only");
  assert.equal(compareWithIce(ko("g"), [{ ok: true, donorName: "", summary: { rsq: 0.9 }, contribs: [] }]).status, "ice_failed", "no score from ICE");
});

test("ICE jobs: one per sample for a knockout, one per donor for a point edit; unsupported samples are listed with a reason", () => {
  const guide = (spacer) => ({ spacer }); const trace = (n) => `/data/${n}.ab1`;
  const tasks = [
    { name: "k1", workflow: "deletion", spec: { guides: [guide("AGGTGCCAGGCAGGACATCG"), guide("AAGGCCACGGGGTGGGTGCG")] }, controls: [trace("wt")], edited: [trace("k1")] },
    { name: "e1", workflow: "snp", spec: { guides: [guide("CGATGACCTGCAGAAGCGCC")], donors: [{ name: "A", sequence: "ACGT".repeat(30) }, { name: "B", sequence: "TTGCA".repeat(20) }] }, controls: [trace("wt")], edited: [trace("e1")] },
    { name: "e2", workflow: "snp", spec: { guides: [guide("CGATGACCTGCAGAAGCGCC")], donors: [{ name: "big", sequence: "ACGT".repeat(100) }] }, controls: [trace("wt")], edited: [trace("e2")] },
    { name: "b1", workflow: "base_edit", spec: { guides: [guide("CGATGACCTGCAGAAGCGCC")] }, controls: [trace("wt")], edited: [trace("b1")] },
    { name: "g1", workflow: "knockout", spec: { guides: [guide("ACGT")] }, controls: [trace("wt")], edited: [trace("g1")] },
  ];
  const { jobs, skipped } = iceJobsFor(tasks);
  assert.deepEqual(jobs.map((j) => [j.sample, j.donorName]), [["k1", ""], ["e1", "A"], ["e1", "B"]]);
  assert.deepEqual(jobs[0].guides, ["AGGTGCCAGGCAGGACATCG", "AAGGCCACGGGGTGGGTGCG"]); assert.equal(jobs[0].donor, ""); assert.equal(jobs[1].donor.length, 120);
  assert.deepEqual(skipped.map((s) => s.sample), ["e2", "b1", "g1"]); assert.match(skipped[0].reason, /20 to 300 bases/); assert.match(skipped[1].reason, /base_edit/); assert.match(skipped[2].reason, /17 to 23/);
  assert.equal(iceJobsFor([{ ...tasks[0], nuclease: "SaCas9" }]).skipped[0].reason, "ICE supports SpCas9 only");
});

test("a run: rows, summary and CSV, with the alleles named", () => {
  const results = [ko("a", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] }), ko("het"), { ok: false, name: "bad" }];
  const byIndex = new Map([[0, [run(99, 0.99, [[-78, 0.99]])]], [1, [run(58, 0.63, [[-78, 0.41]])]]]);
  const rows = compareRun(results, byIndex);
  assert.deepEqual(rows.map((r) => [r.sample, r.status]), [["a", "agree"], ["het", "ice_fit_poor"]]);
  assert.deepEqual(iceSummary(rows), { total: 2, agree: 1, differ: 0, ice_fit_poor: 1, ice_failed: 0, not_run: 0 });
  const lines = iceCsv(rows).trim().split("\n"); assert.equal(lines.length, 3); assert.match(lines[0], /^sample,workflow,assured_qc_score,assured_qc_interval,ice_score/);
  assert.match(lines[2], /het,deletion,100,97-100,58,.*,0\.63,42,ice_fit_poor/); assert.match(lines[2], /-77 \(50%\)/);
});

test("the plot and the report section are valid; the plot converts to a PDF", async () => {
  const rows = compareRun([ko("a", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] }), ko("het"), ko("diff", { contributions: [{ kind: "deletion_between_cuts", size: -78, fraction: 1 }] })],
    new Map([[0, [run(100, 1, [[-78, 1]])]], [1, [run(58, 0.63, [[-78, 0.41]])]], [2, [run(40, 0.95, [[-78, 0.4]])]]]));
  const svg = iceScatterSvg(rows); assert.match(svg, /^<svg /); assert.ok(svg.trimEnd().endsWith("</svg>")); assert.match(svg, />het \(R\u00b2 0\.63\)</); assert.match(svg, />diff \(R\u00b2 0\.95\)</); assert.ok(!/>a</.test(svg.replace(/<title>.*?<\/title>/g, "")), "agreeing clones are counted, not labelled");
  assert.match(svg, /agree, 1 /); assert.match(svg, /ICE fit poor, 1 /); assert.match(svg, /differ, 1/);
  const pdf = await convertSvg(svg, { compress: false }); assert.deepEqual(pdf.unsupported, []); assert.equal(pdf.lost, 0);
  const html = iceSectionHtml(rows, { version: "1.2.0", python: "3.11" }); assert.match(html, /Comparison with ICE/); assert.match(html, /1\.2\.0/); assert.match(html, /not comparable/); assert.equal((html.match(/<tr>/g) || []).length, 4);
  const truth = loadTruth(); const k = truth.apoe_snp_clone_hom; const real = analyseSample({ name: "clone_hom", spec: loadSpec("apoe_r176c_snp"), control: loadTrace(k.control), edited: loadTrace(k.edited), sampleType: "clone" });
  const report = buildQcReportHtml({ workflow: "snp", results: withRecommendations([{ ...real, workflow: "snp", decision: { decision: "accept", reasons: [], next: "" } }]), extraHtml: html }); assert.ok(report.indexOf("Comparison with ICE") > report.indexOf("Summary") && report.indexOf("Comparison with ICE") < report.indexOf("<h2>Samples</h2>"));
  assert.equal(iceScatterSvg([]), ""); assert.equal(iceSectionHtml([]), "");
});

test("ICE results are reused only for the same samples, guides, donors and traces", () => {
  const jobs = [{ sample: "a", donorName: "", guides: ["AAAA"], donor: "", control: "/d/wt.ab1", edited: "/d/a.ab1" }, { sample: "e", donorName: "A", guides: ["AAAA"], donor: "ACGT", control: "/d/wt.ab1", edited: "/d/e.ab1" }];
  const shaOf = (p) => `sha-${path.basename(p)}`; const results = new Map([[0, { ok: true, summary: { ko_score: 99 } }], [1, { ok: true, summary: { hdr_pct: 80 } }]]);
  const records = iceRawRecords(jobs, results, shaOf); assert.equal(records[0].controlSha256, "sha-wt.ab1"); assert.equal(records[1].donorName, "A"); assert.equal(records[0].control, "wt.ab1");
  const ok = matchIceRaw(jobs, records, shaOf); assert.deepEqual(ok.problems, []); assert.equal(ok.results.get(1).summary.hdr_pct, 80);
  assert.match(matchIceRaw(jobs, records.slice(1), shaOf).problems[0], /no earlier ICE result for a/);
  assert.match(matchIceRaw([{ ...jobs[0], guides: ["CCCC"] }, jobs[1]], records, shaOf).problems[0], /other guides or another donor/);
  assert.match(matchIceRaw(jobs, records, (p) => `x-${path.basename(p)}`).problems[0], /other traces/);
});

test("JSON written by Python (NaN, Infinity) is read; the words inside strings are left alone", () => {
  const text = '{"a": NaN, "b": [1, -Infinity, 2], "c": "NaN in a string, Infinity", "d": {"x": Infinity, "y": 0.5}}';
  assert.deepEqual(parseIceJson(text), { a: null, b: [1, null, 2], c: "NaN in a string, Infinity", d: { x: null, y: 0.5 } });
});

test("ICE check: a missing checkout, a Python that cannot start, one without the packages, and a good one", () => {
  const dir = mkdtempSync(path.join(tmpdir(), "ice-")); assert.match(checkIce(dir, "python").error, /not a Synthego ICE checkout/);
  mkdirSync(path.join(dir, "ice")); writeFileSync(path.join(dir, "ice", "analysis.py"), ""); writeFileSync(path.join(dir, "ice", "__version__.py"), "__version__ = '1.2.0'\n");
  assert.match(checkIce(dir, "python", () => ({ error: { code: "EPERM" } })).error, /could not start python \(EPERM\)/);
  assert.match(checkIce(dir, "python", () => ({ status: 1, stderr: "Traceback\nModuleNotFoundError: No module named 'Bio'\n" })).error, /No module named 'Bio'/);
  assert.deepEqual(checkIce(dir, "python", () => ({ status: 0, stdout: "3.11.9\n" })), { ok: true, version: "1.2.0", python: "3.11.9" });
});

// a stand-in for the Python process: reads the jobs file it is given and writes what run_ice_jobs.py would
const fakeSpawn = (behaviour) => (cmd, args) => {
  const handlers = {}; const child = { stdout: { on() {} }, stderr: { on: (e, f) => { child.stderrFn = f; } }, on: (e, f) => { handlers[e] = f; } };
  setImmediate(() => {
    const [script, iceDir, jobsFile, outFile] = args; assert.equal(script, ICE_SCRIPT); const jobs = JSON.parse(readFileSync(jobsFile, "utf8"));
    if (behaviour === "ok") { const out = {}; jobs.forEach((j, k) => { out[k] = { ok: true, summary: { ko_score: 90 + k, rsq: 0.9, mean_discord_after: "@NAN@" }, contribs: [{ rel: 0.5, indel: -78 }] }; }); writeFileSync(outFile, JSON.stringify(out).replaceAll('"@NAN@"', "NaN")); handlers.close(0); }
    else if (behaviour === "crash") { child.stderrFn("Traceback\nImportError: boom\n"); handlers.close(1); }
    else handlers.error(new Error("spawn EPERM"));
  });
  return child;
};
test("the runner shares jobs over processes, reads their results, and reports a process that fails", async () => {
  const jobs = Array.from({ length: 5 }, (_, i) => ({ control: `/d/wt.ab1`, edited: `/d/c${i}.ab1`, guides: ["AAAA"], donor: "" })); const seen = [];
  const dir = path.join(mkdtempSync(path.join(tmpdir(), "ice-run-")), "w");
  const ok = await runIceJobs({ iceDir: "/ice", python: "python", jobs, workdir: dir, concurrency: 2, spawnFn: fakeSpawn("ok"), onProgress: (d, n) => seen.push([d, n]) });
  assert.equal(ok.size, 5); assert.ok([...ok.values()].every((r) => r.ok && r.summary.mean_discord_after === null)); assert.equal(ok.get(0).summary.ko_score, 90); assert.equal(ok.get(2).summary.ko_score, 91, "jobs 0, 2, 4 went to the first process"); assert.deepEqual(seen[seen.length - 1], [5, 5]);
  const crash = await runIceJobs({ iceDir: "/ice", python: "python", jobs: jobs.slice(0, 2), workdir: path.join(dir, "c"), concurrency: 1, spawnFn: fakeSpawn("crash") });
  assert.ok([...crash.values()].every((r) => !r.ok && /ImportError: boom/.test(r.error)));
  const blocked = await runIceJobs({ iceDir: "/ice", python: "python", jobs: jobs.slice(0, 1), workdir: path.join(dir, "b"), concurrency: 1, spawnFn: fakeSpawn("blocked") });
  assert.match(blocked.get(0).error, /could not start python: spawn EPERM/);
});
