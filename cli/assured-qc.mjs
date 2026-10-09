#!/usr/bin/env node
// Assured QC command line: a folder of .ab1 files (or a sample sheet) in; results, a report, a plate view and a provenance record out.
//
//   node cli/assured-qc.mjs run <folder> [--sheet sheet.csv] [--design design.json|report.html] [--guides G1,G2] [--donor SEQ]
//        [--workflow auto|knockout|deletion|snp|tag_small|base_edit] [--nuclease SpCas9] [--editor ABE|CBE --window 4-8 --target 6]
//        [--type clone|pool] [--goal homozygous|any] [--threads N] [--out results/]
//   node cli/assured-qc.mjs figure --control control.ab1 --edited edited.ab1 [--design design.json | --guides G1,G2] [--out figure.svg]
//
// Nothing is uploaded anywhere: the files are read from disk and analysed in this process and its worker threads.
import { readdirSync, readFileSync, writeFileSync, mkdirSync, statSync, existsSync } from "node:fs";
import { Worker } from "node:worker_threads";
import { cpus } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { parseAbif } from "../src/abif.js";
import { planDesign, WORKFLOWS } from "../src/batch.js";
import { autoPair, parseSampleSheet, resolveSheet } from "../src/sheet.js";
import { parseDesignReport } from "../src/designReport.js";
import { buildQcReportHtml, buildResultsCsv, QC_VERSION } from "../src/report.js";
import { buildPlateReportHtml, resequenceListCsv } from "../src/plateReport.js";
import { tally, wellOf } from "../src/decision.js";
import { chromatogramFigure } from "../src/figure.js";

const here = path.dirname(fileURLToPath(import.meta.url));

function parseArgs(argv) {
  const args = { _: [] };
  for (let i = 0; i < argv.length; i += 1) { const a = argv[i]; if (a.startsWith("--")) { const key = a.slice(2); const next = argv[i + 1]; if (next === undefined || next.startsWith("--")) args[key] = true; else { args[key] = next; i += 1; } } else args._.push(a); }
  return args;
}
const walk = (dir) => readdirSync(dir).flatMap((name) => { const p = path.join(dir, name); return statSync(p).isDirectory() ? walk(p) : [p]; });
const readTrace = (p) => { const b = readFileSync(p); return parseAbif(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength), path.basename(p)); };

export async function main(argv) {
  const args = parseArgs(argv); const command = args._[0];
  if (!command || command === "help" || args.help) { console.log(readFileSync(fileURLToPath(import.meta.url), "utf8").split("\n").slice(2, 9).map((l) => l.replace(/^\/\/ ?/, "")).join("\n")); return 0; }
  if (command === "figure") {
    if (!args.control || !args.edited) { console.error("figure needs --control FILE --edited FILE (and --guides or --design)."); return 1; }
    const control = readTrace(args.control); const edited = readTrace(args.edited);
    const guides = args.guides ? String(args.guides).split(/[,;\s]+/).filter(Boolean) : [];
    const designSpec = args.design && /^\s*[{[]/.test(readFileSync(args.design, "utf8")) ? JSON.parse(readFileSync(args.design, "utf8")) : null;
    const plan = planDesign({ designSpec, control, guides, donor: args.donor || "", gene: "sample", nuclease: args.nuclease || "SpCas9", editedForCut: [edited] });
    if (plan.error) { console.error(plan.error); return 1; }
    const { svg, result } = chromatogramFigure({ name: args.name || path.basename(args.edited).replace(/\.ab1$/i, ""), controlTrace: control, editedTrace: edited, spec: plan.spec, sampleType: args.type || "clone", workflow: args.workflow || (plan.spec.donors?.length ? "snp" : plan.spec.guides.length > 1 ? "deletion" : "knockout") });
    const target = args.out || "figure.svg"; writeFileSync(target, svg); console.log(`${target}${result && result.decision ? `: ${result.decision.decision}` : ""}`); return 0;
  }
  if (command !== "run") { console.error(`Unknown command ${command}. Use: run <folder> | figure --control FILE --edited FILE`); return 1; }
  const folder = args._[1]; if (!folder || !existsSync(folder)) { console.error("Give the folder with the .ab1 files."); return 1; }
  const out = args.out || path.join(folder, "assured-qc-results"); mkdirSync(out, { recursive: true });
  const files = walk(folder).filter((p) => /\.(ab1|abi)$/i.test(p) && !p.startsWith(path.resolve(out)));
  const byName = new Map(files.map((p) => [path.basename(p), p]));
  // 1. pairing: a sample sheet if given, otherwise by names
  let jobs = []; const problems = [];
  if (args.sheet) {
    const sheet = parseSampleSheet(readFileSync(args.sheet, "utf8")); problems.push(...sheet.errors);
    const resolved = resolveSheet(sheet.rows, [...byName.keys()]); problems.push(...resolved.problems);
    jobs = resolved.jobs.filter((j) => j.ok).map((j) => ({ ...j, confidence: "sheet" }));
  } else {
    const paired = autoPair([...byName.keys()]);
    jobs = paired.pairs.map((p) => ({ sample: p.sample, edited: p.edited, control: p.control, type: p.type, confidence: p.confidence, reason: p.reason, guides: [], donor: "", workflow: "", nuclease: "" }));
    paired.unpaired.forEach((u) => problems.push(`${u.name}: not analysed (${u.reason})`)); paired.problems.forEach((x) => problems.push(`Pairing to check: ${x}`));
  }
  if (!jobs.length) { console.error(`No samples to analyse.\n${problems.join("\n")}`); return 1; }
  // 2. design: file, guides on the command line, or (per sheet row) guides in the sheet; else inferred from the traces
  let designSpec = null; let design = { guides: args.guides ? String(args.guides).split(/[,;\s]+/).filter(Boolean) : [], donor: args.donor || "" };
  if (args.design) {
    const text = readFileSync(args.design, "utf8");
    if (/^\s*[{[]/.test(text)) designSpec = JSON.parse(text);
    else { const designs = parseDesignReport(text); const chosen = designs[Number(args["design-index"] || 0)]; design = { guides: chosen.guides.map((g) => g.sequence), donor: (chosen.donors.find((d) => d.recommended) || chosen.donors[0] || {}).sequence || "" }; console.log(`Design report: ${chosen.title}; ${chosen.guides.length} guide(s), ${chosen.donors.length} donor(s).`); }
  }
  // 3. one spec per distinct control + guides + donor (the trace is only read once per group)
  const cache = new Map(); const tasks = []; const nuclease = args.nuclease || "SpCas9";
  const baseEditor = args.editor ? { editor: args.editor, window: String(args.window || "4-8").split("-").map(Number), target: args.target ? Number(args.target) : null } : null;
  jobs.forEach((job, index) => {
    const guides = job.guides?.length ? job.guides : design.guides; const donor = job.donor || design.donor;
    const controlPath = byName.get(job.control[0]); const key = JSON.stringify([controlPath, guides, donor, job.nuclease || nuclease]);
    if (!cache.has(key)) {
      const control = readTrace(controlPath);
      const plan = planDesign({ designSpec, control, guides, donor, gene: path.basename(folder), nuclease: job.nuclease || nuclease, baseEditor, editedForCut: jobs.filter((j) => j.control[0] === job.control[0]).slice(0, 6).map((j) => { try { return readTrace(byName.get(j.edited[0])); } catch { return null; } }).filter(Boolean) });
      cache.set(key, plan);
    }
    const plan = cache.get(key);
    const workflow = args.workflow && args.workflow !== "auto" ? args.workflow : baseEditor ? "base_edit" : job.workflow || (plan.spec && (plan.spec.donors || []).length ? ((plan.spec.donors[0].insertBp || 0) > 0 ? "tag_small" : "snp") : (plan.spec && plan.spec.guides.length > 1 ? "deletion" : "knockout"));
    if (plan.error) { problems.push(`${job.sample}: ${plan.error}`); return; }
    tasks.push({ index, name: job.sample, spec: plan.spec, workflow, sampleType: args.type || job.type || "clone", controls: job.control.map((n) => byName.get(n)), edited: job.edited.map((n) => byName.get(n)), options: { library: {}, baseEdit: {} }, decisionOptions: { goal: args.goal === "any" ? "any" : "homozygous" }, meta: { confidence: job.confidence, reason: job.reason, planNotes: plan.warnings || [] } });
  });
  // 4. analysis on a pool of worker threads
  const threads = Math.max(1, Math.min(Number(args.threads) || Math.max(1, cpus().length - 1), tasks.length)); const started = Date.now();
  const results = new Array(tasks.length); let next = 0; let done = 0;
  await new Promise((resolve) => {
    const finish = () => { if (done === tasks.length) { workers.forEach((w) => w.terminate()); resolve(); } };
    const workers = Array.from({ length: threads }, () => { const w = new Worker(path.join(here, "worker.mjs")); w.on("message", ({ index, result }) => { results[index] = result; done += 1; if (done % 20 === 0 || done === tasks.length) process.stderr.write(`\r${done}/${tasks.length} samples`); if (next < tasks.length) { w.postMessage(tasks[next]); next += 1; } finish(); }); return w; });
    // dispatch: each task carries its own index into `results` (its position in tasks)
    tasks.forEach((t, i) => { t.index = i; });
    workers.forEach((w) => { if (next < tasks.length) { w.postMessage(tasks[next]); next += 1; } });
  });
  process.stderr.write("\n");
  results.forEach((r, i) => { r.name = tasks[i].name; r.pairing = tasks[i].meta; r.well = wellOf(path.basename(tasks[i].edited[0])); r.workflow = tasks[i].workflow; if (tasks[i].meta.planNotes.length) r.warnings = [...tasks[i].meta.planNotes, ...(r.warnings || [])]; });
  const seconds = (Date.now() - started) / 1000;
  // 5. outputs
  const workflowKey = results[0].workflow; const csv = buildResultsCsv(workflowKey, results);
  writeFileSync(path.join(out, "results.csv"), csv);
  const slim = results.map(({ charts, ...rest }) => rest); writeFileSync(path.join(out, "results.json"), JSON.stringify(slim, null, 1));
  writeFileSync(path.join(out, "report.html"), buildQcReportHtml({ title: `Sanger QC report: ${path.basename(path.resolve(folder))}`, workflow: workflowKey, design: design.guides.join(", "), results, notes: problems }));
  writeFileSync(path.join(out, "plate.html"), buildPlateReportHtml({ title: path.basename(path.resolve(folder)), results }));
  writeFileSync(path.join(out, "resequence.csv"), resequenceListCsv(results));
  writeFileSync(path.join(out, "run.json"), JSON.stringify({ tool: "Assured QC", version: QC_VERSION, finished: new Date().toISOString(), seconds, threads, samples: results.length, samplesPerHour: Math.round((results.length / seconds) * 3600), arguments: { ...args, _: undefined }, workflow: workflowKey, nuclease, pairingProblems: problems, decisions: tally(results.map((r) => r.decision)) }, null, 1));
  const t = tally(results.map((r) => r.decision)); console.log(`${results.length} samples in ${seconds.toFixed(1)} s on ${threads} thread(s) (${Math.round((results.length / seconds) * 3600)} per hour): accept ${t.accept}, hold ${t.hold}, re-sequence ${t["re-sequence"]}, reject ${t.reject}. Results in ${out}`);
  problems.forEach((p) => console.log(`  note: ${p}`));
  return 0;
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) main(process.argv.slice(2)).then((code) => process.exit(code));
