#!/usr/bin/env node
// Assured QC command line: a folder of .ab1 files (or a sample sheet) in; results, a report, a plate view and a provenance record out.
//
//   node cli/assured-qc.mjs run <folder> [--sheet sheet.csv] [--design design.json|report.html] [--guides G1,G2] [--donor SEQ]
//        [--workflow auto|knockout|deletion|snp|tag_small|base_edit] [--nuclease SpCas9] [--editor ABE|CBE --window 4-8 --target 6]
//        [--type clone|pool] [--goal homozygous|any] [--threads N] [--out results/]
//        [--genbank gene.gb [--transcript NM_000000] [--regions "TIR domain:172-314,Death domain:67-122"] [--epitope 200-250]]   (functional knockout assessment)
//   node cli/assured-qc.mjs figure --control control.ab1 --edited edited.ab1 [--design design.json | --guides G1,G2] [--out figure.svg|figure.pdf]
//   node cli/assured-qc.mjs recommend results.json [--observations protein.csv] [--panel N] [--out recommendation.csv]   (rank the clones, choose the panel to confirm)
//   node cli/assured-qc.mjs calibrate results.json observations.csv [--out calibration.md]   (how often the functional call was borne out)
//   run also takes --observations protein.csv (sample, protein[, rna]), --panel N (clones to carry forward first: 3 for a knockout, 6 for an edit) --clone-figures (one annotated chromatogram per clone, aligned with the control) and --ice FOLDER (also run Synthego ICE, a checkout of github.com/synthego-open/ice, and compare: ice_comparison.csv, ice_vs_assuredqc.svg, ice_raw.json and a section in the report; --ice-python PYTHON, --ice-jobs N; --ice-raw ice_raw.json reuses earlier results for the same traces; --ice-plan writes the ICE jobs to ice_jobs.json for scripts/run_ice_plan.py; --ice-results ICE_export.zip compares with results downloaded from the ICE web tool (the zip, its summary .csv or the unzipped folder) without running ICE; ICE changes no call or rank) and --no-pdf (every figure is also written as a vector PDF next to its SVG; put this option last)
//
// Nothing is uploaded anywhere: the files are read from disk and analysed in this process and its worker threads.
import { readdirSync, readFileSync, writeFileSync, mkdirSync, statSync, existsSync } from "node:fs";
import { Worker } from "node:worker_threads";
import { cpus } from "node:os";
import { createHash } from "node:crypto";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { parseAbif } from "../src/abif.js";
import { planDesign, WORKFLOWS, headline } from "../src/batch.js";
import { parseCodingFields } from "../src/coding.js";
import { calibrate, parseObservations, calibrationMarkdown } from "../src/calibration.js";
import { withRecommendations, applyObservations, recommendationCsv, panelOf } from "../src/recommend.js";
import { runProteinFigure } from "../src/proteinFigure.js";
import { runSummaryFigure } from "../src/summaryFigure.js";
import { autoPair, parseSampleSheet, resolveSheet } from "../src/sheet.js";
import { parseDesignReport } from "../src/designReport.js";
import { reportDonors } from "../src/quick.js";
import { cloneSelectionFigure } from "../src/selectionFigure.js";
import { svgToPdf } from "../src/svgPdf.js";
import { iceJobsFor, compareRun, iceCsv, iceSectionHtml, iceScatterSvg, iceSummary, iceRawRecords, matchIceRaw, parseIceJson, ICE_FIT_OK, ICE_DIFFERENCE } from "../src/iceCompare.js";
import { checkIce, icePython, runIceJobs } from "./ice-runner.mjs";
import { readIceExport, compareWithExport } from "../src/iceImport.js";
import { buildQcReportHtml, buildResultsCsv, QC_VERSION } from "../src/report.js";
import { buildPlateReportHtml, resequenceListCsv } from "../src/plateReport.js";
import { tally, wellOf } from "../src/decision.js";
import { chromatogramFigure } from "../src/figure.js";
import { alignedGroups, cloneAlignedFigure, clonesAlignedFigure, clonesAlignmentFigure } from "../src/alignedFigure.js";

const here = path.dirname(fileURLToPath(import.meta.url));

// --donor takes one sequence, or several separated by commas (used together, as the ssODNs of one design)
const donorArg = (value) => { const parts = String(value || "").split(/[,;\s]+/).filter(Boolean); return parts.length > 1 ? parts : (parts[0] || ""); };

function parseArgs(argv) {
  const args = { _: [] };
  for (let i = 0; i < argv.length; i += 1) { const a = argv[i]; if (a.startsWith("--")) { const key = a.slice(2); const next = argv[i + 1]; if (next === undefined || next.startsWith("--")) args[key] = true; else { args[key] = next; i += 1; } } else args._.push(a); }
  return args;
}
const walk = (dir) => readdirSync(dir).flatMap((name) => { const p = path.join(dir, name); return statSync(p).isDirectory() ? walk(p) : [p]; });
const readTrace = (p) => { const b = readFileSync(p); return parseAbif(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength), path.basename(p)); };

export async function main(argv) {
  const args = parseArgs(argv); const command = args._[0];
  if (!command || command === "help" || args.help) { console.log(readFileSync(fileURLToPath(import.meta.url), "utf8").split("\n").slice(2, 14).map((l) => l.replace(/^\/\/ ?/, "")).join("\n")); return 0; }
  if (command === "calibrate") {
    // assured-qc calibrate results.json observations.csv [--out calibration.md]: how often the functional knockout call was borne out by protein (and RNA) results
    const [resultsPath, observationsPath] = [args._[1], args._[2]];
    if (!resultsPath || !observationsPath || !existsSync(resultsPath) || !existsSync(observationsPath)) { console.error("calibrate needs results.json (from a run with a gene record) and a table of protein results (sample, protein[, rna])."); return 1; }
    const parsed = parseObservations(readFileSync(observationsPath, "utf8")); if (parsed.errors.length) { parsed.errors.forEach((e) => console.error(e)); return 1; }
    const result = calibrate(JSON.parse(readFileSync(resultsPath, "utf8")), parsed.rows); const markdown = calibrationMarkdown(result);
    const target = args.out || "calibration.md"; writeFileSync(target, markdown); writeFileSync(target.replace(/\.md$/i, "") + ".json", JSON.stringify(result, null, 1));
    console.log(`${result.n} clones scored; predicted knockouts ${result.knockoutCall.called}, protein absent in ${result.knockoutCall.confirmed}. ${result.note} Written to ${target}.`); return 0;
  }
  if (command === "recommend") {
    // assured-qc recommend results.json [--observations protein.csv] [--panel N] [--out recommendation.csv]: rank the clones and name the panel to take to the bench
    const resultsPath = args._[1]; if (!resultsPath || !existsSync(resultsPath)) { console.error("recommend needs results.json (from a run)."); return 1; }
    let results = JSON.parse(readFileSync(resultsPath, "utf8")); const notes = [];
    if (args.observations) {
      if (!existsSync(args.observations)) { console.error(`Protein results not found: ${args.observations}`); return 1; }
      const parsed = parseObservations(readFileSync(args.observations, "utf8")); if (parsed.errors.length) { parsed.errors.forEach((e) => console.error(e)); return 1; }
      const applied = applyObservations(results, parsed.rows); results = applied.results; notes.push(...applied.problems); console.log(`Protein results applied to ${applied.matched} sample(s).`);
    }
    results = withRecommendations(results, { panel: Number(args.panel) || undefined });
    const target = args.out || "recommendation.csv"; writeFileSync(target, recommendationCsv(results)); writeFileSync(target.replace(/\.csv$/i, "") + ".json", JSON.stringify(results.map((r) => ({ name: r.name, decision: r.decision, confirmation: r.confirmation || null, recommendation: r.recommendation })), null, 1));
    console.log(`Take forward first: ${panelOf(results).join(", ") || "none"}. Written to ${target}.`); notes.forEach((n) => console.log(`  note: ${n}`)); return 0;
  }
  if (command === "figure") {
    if (!args.control || !args.edited) { console.error("figure needs --control FILE --edited FILE (and --guides or --design)."); return 1; }
    const control = readTrace(args.control); const edited = readTrace(args.edited);
    const guides = args.guides ? String(args.guides).split(/[,;\s]+/).filter(Boolean) : [];
    const designSpec = args.design && /^\s*[{[]/.test(readFileSync(args.design, "utf8")) ? JSON.parse(readFileSync(args.design, "utf8")) : null;
    const plan = planDesign({ designSpec, control, guides, donor: donorArg(args.donor), gene: "sample", nuclease: args.nuclease || "SpCas9", editedForCut: [edited] });
    if (plan.error) { console.error(plan.error); return 1; }
    const { svg, result } = chromatogramFigure({ name: args.name || path.basename(args.edited).replace(/\.ab1$/i, ""), controlTrace: control, editedTrace: edited, spec: plan.spec, sampleType: args.type || "clone", workflow: args.workflow || (plan.spec.donors?.length ? "snp" : plan.spec.guides.length > 1 ? "deletion" : "knockout") });
    const target = args.out || "figure.svg"; writeFileSync(target, /\.pdf$/i.test(target) ? await svgToPdf(svg) : svg); console.log(`${target}${result && result.decision ? `: ${result.decision.decision}` : ""}`); return 0;
  }
  if (command !== "run") { console.error(`Unknown command ${command}. Use: run <folder> | figure --control FILE --edited FILE | calibrate results.json observations.csv | recommend results.json`); return 1; }
  const folder = args._[1]; if (!folder || !existsSync(folder)) { console.error("Give the folder with the .ab1 files."); return 1; }
  const out = args.out || path.join(folder, "assured-qc-results"); mkdirSync(out, { recursive: true });
  // optional comparison with Synthego ICE (a source checkout given with --ice); checked now so that a wrong folder or Python is reported before the analysis
  // --ice-raw: reuse the ice_raw.json of an earlier run (same samples, guides, donors and traces) instead of running ICE again
  const iceRawPath = args["ice-raw"]; if (iceRawPath && (iceRawPath === true || !existsSync(iceRawPath))) { console.error(`ICE results not found: ${iceRawPath}`); return 1; }
  // --ice-results: results downloaded from the ICE web tool (its zip, the summary .csv, or the unzipped folder); read now so that a wrong file is reported before the analysis
  const iceResultsPath = args["ice-results"]; let iceExport = null; let iceExportFiles = [];
  if (iceResultsPath) {
    if (iceResultsPath === true || !existsSync(iceResultsPath)) { console.error(`ICE results not found: ${iceResultsPath}`); return 1; }
    if (args.ice || args["ice-raw"] || args["ice-plan"]) { console.error("--ice-results (results downloaded from the ICE web tool) cannot be combined with --ice, --ice-raw or --ice-plan."); return 1; }
    const given = statSync(iceResultsPath).isDirectory() ? readdirSync(iceResultsPath).filter((n) => /^summary.*\.csv$/i.test(n) || (/\.zip$/i.test(n) && !readdirSync(iceResultsPath).some((m) => /^summary.*\.csv$/i.test(m)))).map((n) => path.join(iceResultsPath, n)) : [iceResultsPath];
    if (!given.length) { console.error(`ICE results: no summary .csv or .zip in ${iceResultsPath}`); return 1; }
    iceExportFiles = given.map((p) => ({ name: path.basename(p), bytes: new Uint8Array(readFileSync(p)) }));
    iceExport = await readIceExport(iceExportFiles);
    if (!iceExport.rows.length) { console.error(`ICE results: ${iceExport.problems.join(" ") || "no samples found"}`); return 1; }
    iceExport.problems.forEach((p) => console.error(`ICE results: ${p}`));
  }
  const python = icePython(args); const iceCheck = args.ice && !iceRawPath ? checkIce(args.ice, python) : null;
  if (iceCheck && !iceCheck.ok) { console.error(iceCheck.error); return 1; }
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
  let designSpec = null; let design = { guides: args.guides ? String(args.guides).split(/[,;\s]+/).filter(Boolean) : [], donor: donorArg(args.donor) };
  if (args.design) {
    const text = readFileSync(args.design, "utf8");
    if (/^\s*[{[]/.test(text)) designSpec = JSON.parse(text);
    else { const designs = parseDesignReport(text); const chosen = designs[Number(args["design-index"] || 0)]; design = { guides: chosen.guides.map((g) => g.sequence), donor: reportDonors(chosen.donors) }; console.log(`Design report: ${chosen.title}; ${chosen.guides.length} guide(s), ${chosen.donors.length} donor(s)${Array.isArray(design.donor) ? `; the ${design.donor.length} recommended donors are used together` : ""}.`); }
  }
  // 3. one spec per distinct control + guides + donor (the trace is only read once per group)
  // the gene record, if given: each allele is translated and the knockout is judged on the protein (NMD, truncation, splice sites, in-frame changes)
  let genbank = "";
  const parsedCoding = parseCodingFields({ regionsText: args.regions || "", epitopeText: args.epitope || "", gene: args.gene || "", transcript: args.transcript || "" });
  if (parsedCoding.problems.length) { parsedCoding.problems.forEach((p) => console.error(p)); return 1; }
  const codingOptions = parsedCoding.options;
  if (args.genbank) {
    if (!existsSync(args.genbank)) { console.error(`Gene record not found: ${args.genbank}`); return 1; }
    genbank = readFileSync(args.genbank, "utf8");
  }
  const cache = new Map(); const tasks = []; const nuclease = args.nuclease || "SpCas9";
  const baseEditor = args.editor ? { editor: args.editor, window: String(args.window || "4-8").split("-").map(Number), target: args.target ? Number(args.target) : null } : null;
  jobs.forEach((job, index) => {
    const guides = job.guides?.length ? job.guides : design.guides; const donor = job.donor || design.donor;
    const controlPath = byName.get(job.control[0]); const key = JSON.stringify([controlPath, guides, donor, job.nuclease || nuclease]);
    if (!cache.has(key)) {
      const control = readTrace(controlPath);
      const plan = planDesign({ designSpec, control, guides, donor, gene: path.basename(folder), nuclease: job.nuclease || nuclease, baseEditor, genbank, codingOptions, editedForCut: jobs.filter((j) => j.control[0] === job.control[0]).slice(0, 6).map((j) => { try { return readTrace(byName.get(j.edited[0])); } catch { return null; } }).filter(Boolean) });
      cache.set(key, plan);
    }
    const plan = cache.get(key);
    const workflow = args.workflow && args.workflow !== "auto" ? args.workflow : baseEditor ? "base_edit" : job.workflow || (plan.spec && (plan.spec.donors || []).length ? ((plan.spec.donors[0].insertBp || 0) > 0 ? "tag_small" : "snp") : (plan.spec && plan.spec.guides.length > 1 ? "deletion" : "knockout"));
    if (plan.error) { problems.push(`${job.sample}: ${plan.error}`); return; }
    tasks.push({ index, name: job.sample, spec: plan.spec, workflow, sampleType: args.type || job.type || "clone", controls: job.control.map((n) => byName.get(n)), edited: job.edited.map((n) => byName.get(n)), options: { library: {}, baseEdit: {}, codingInput: plan.codingInput || null }, decisionOptions: { goal: args.goal === "any" ? "any" : "homozygous" }, meta: { confidence: job.confidence, reason: job.reason, planNotes: plan.warnings || [] } });
  });
  if ((args.regions || args.epitope) && !tasks.some((t) => t.options.codingInput)) problems.push("--regions or --epitope were given but no gene record is available (use --genbank, or a design file with a coding block); they were ignored.");
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
  // protein results from the bench change the decision of the samples they cover; every clone is then ranked and the panel to confirm is named
  if (args.observations) {
    if (!existsSync(args.observations)) { console.error(`Protein results not found: ${args.observations}`); return 1; }
    const parsed = parseObservations(readFileSync(args.observations, "utf8")); if (parsed.errors.length) { parsed.errors.forEach((e) => console.error(e)); return 1; }
    const applied = applyObservations(results, parsed.rows); applied.results.forEach((r, i) => { results[i] = r; }); applied.problems.forEach((p) => problems.push(`Protein results: ${p}`));
  }
  withRecommendations(results, { panel: Number(args.panel) || undefined }).forEach((r, i) => { results[i] = r; });
  // 4b. comparison with ICE: run on the same traces (one process per share of the jobs) or taken from an earlier run; it changes no call, rank or score above
  let iceRows = null; let iceMeta = null; const icePlan = Boolean(args["ice-plan"]) && !iceCheck && !iceRawPath;
  if (iceCheck || iceRawPath || icePlan) {
    const { jobs: iceJobs, skipped } = iceJobsFor(tasks, { nuclease }); const iceStarted = Date.now(); const sha = (p) => createHash("sha256").update(readFileSync(p)).digest("hex");
    let ran; let iceVersion = iceCheck ? iceCheck.version : ""; let icePythonVersion = iceCheck ? iceCheck.python : "";
    if (icePlan) {
      // the jobs only: ICE is run elsewhere (scripts/run_ice_plan.py) and its ice_raw.json is given back with --ice-raw
      writeFileSync(path.join(out, "ice_jobs.json"), JSON.stringify({ tool: "Synthego ICE", jobs: iceJobs.map(({ index, workflow, ...job }) => ({ ...job, control: path.resolve(job.control), edited: path.resolve(job.edited) })), notRun: skipped }, null, 1));
      console.log(`ICE jobs (${iceJobs.length}) written to ice_jobs.json. Run: python scripts/run_ice_plan.py <ICE folder> ${path.join(out, "ice_jobs.json")} ${path.join(out, "ice_raw.json")}  then add --ice-raw ${path.join(out, "ice_raw.json")}`);
    } else if (iceRawPath) {
      const raw = parseIceJson(readFileSync(iceRawPath, "utf8")); const matched = matchIceRaw(iceJobs, raw.records, sha);
      if (matched.problems.length) { matched.problems.forEach((p) => console.error(`ICE results: ${p}`)); return 1; }
      ran = matched.results; iceVersion = raw.version || ""; icePythonVersion = raw.python || "";
      console.log(`ICE ${iceVersion}: results of ${iceJobs.length} run(s) reused from ${path.basename(iceRawPath)}.`);
    } else {
      console.log(`ICE ${iceVersion}: ${iceJobs.length} run(s) for ${new Set(iceJobs.map((j) => j.index)).size} of ${tasks.length} sample(s)${skipped.length ? `; ${skipped.length} not run (${[...new Set(skipped.map((s) => s.reason))].join("; ")})` : ""}.`);
      ran = await runIceJobs({ iceDir: args.ice, python, jobs: iceJobs, workdir: path.join(out, "ice_work"), concurrency: Number(args["ice-jobs"]) || threads, onProgress: (d, n) => process.stderr.write(`\rICE ${d}/${n}`) });
      process.stderr.write("\n");
    }
    if (!icePlan) {
    const byIndex = new Map(); iceJobs.forEach((j, k) => { if (!byIndex.has(j.index)) byIndex.set(j.index, []); byIndex.get(j.index).push({ donorName: j.donorName, ...ran.get(k) }); });
    iceRows = compareRun(results, byIndex, skipped);
    writeFileSync(path.join(out, "ice_raw.json"), JSON.stringify({ tool: "Synthego ICE", version: iceVersion, python: icePythonVersion, records: iceRawRecords(iceJobs, ran, sha) }, null, 1));
    iceMeta = { tool: "Synthego ICE", version: iceVersion, source: iceRawPath ? `results reused from ${path.basename(iceRawPath)}` : path.resolve(args.ice), python: icePythonVersion, pythonExecutable: iceRawPath ? "" : python, seconds: Math.round((Date.now() - iceStarted) / 100) / 10, runs: iceJobs.length, nuclease: "SpCas9", fitThreshold: ICE_FIT_OK, differenceThreshold: ICE_DIFFERENCE, donors: "one per ICE run; the best donor is shown", summary: iceSummary(iceRows), rawResults: "ice_raw.json (inputs by SHA-256, as returned by ICE)", note: "ICE is a comparator: it changes no call, rank or score. Used under its research and non-profit licence; not redistributed." };
    }
  }
  if (iceExport) {
    // results from the ICE web tool: matched to the samples by name; nothing is run
    const cmp = compareWithExport({ tasks, results, rows: iceExport.rows }); iceRows = cmp.rows; cmp.problems.forEach((p) => console.error(`ICE results: ${p}`));
    const shaOf = (bytes) => createHash("sha256").update(bytes).digest("hex");
    iceMeta = { tool: "ICE web tool (exported results)", imported: true, source: iceExport.sources.join(", "), files: iceExportFiles.map((f) => ({ file: f.name, sha256: shaOf(f.bytes) })), labelsRead: iceExport.rows.length, labelsNotMatched: cmp.unmatched.map((r) => r.label), fitThreshold: ICE_FIT_OK, differenceThreshold: ICE_DIFFERENCE, donors: "one ICE run per donor for a point edit or tag; the best donor is shown", summary: iceSummary(iceRows), note: "ICE was run by the user in the web tool, not by Assured QC. ICE is a comparator: it changes no call, rank or score." };
    writeFileSync(path.join(out, "ice_results_used.json"), JSON.stringify({ tool: iceMeta.tool, files: iceMeta.files, rows: iceExport.rows }, null, 1));
  }
  // 5. outputs
  const workflowKey = results[0].workflow; const csv = buildResultsCsv(workflowKey, results);
  writeFileSync(path.join(out, "results.csv"), csv);
  const slim = results.map(({ charts, view, ...rest }) => rest); writeFileSync(path.join(out, "results.json"), JSON.stringify(slim, null, 1));
  writeFileSync(path.join(out, "report.html"), buildQcReportHtml({ title: `Sanger QC report: ${path.basename(path.resolve(folder))}`, workflow: workflowKey, design: design.guides.join(", "), results, notes: problems, extraHtml: iceRows ? iceSectionHtml(iceRows, iceMeta) : "" }));
  writeFileSync(path.join(out, "plate.html"), buildPlateReportHtml({ title: path.basename(path.resolve(folder)), results }));
  writeFileSync(path.join(out, "resequence.csv"), resequenceListCsv(results));
  writeFileSync(path.join(out, "recommendation.csv"), recommendationCsv(results));
  // every figure is written as SVG and, unless --no-pdf, as a vector PDF next to it (the PDFs are made after the SVGs, below)
  const pdfJobs = []; const writeFigure = (file, svg) => { writeFileSync(file, svg); if (!args["no-pdf"]) pdfJobs.push([file.replace(/\.svg$/i, ".pdf"), svg]); };
  if (iceRows) { writeFileSync(path.join(out, "ice_comparison.csv"), iceCsv(iceRows)); const iceSvg = iceScatterSvg(iceRows); if (iceSvg) writeFigure(path.join(out, "ice_vs_assuredqc.svg"), iceSvg); }
  const summarySvg = runSummaryFigure(results, { title: `Run summary: ${path.basename(path.resolve(folder))}` }); if (summarySvg) writeFigure(path.join(out, "summary.svg"), summarySvg);
  const productsSvg = runProteinFigure(results); if (productsSvg) writeFigure(path.join(out, "protein_products.svg"), productsSvg);
  const selectionSvg = cloneSelectionFigure(results); if (selectionSvg) writeFigure(path.join(out, "clone_selection.svg"), selectionSvg);
  // all clones laid onto the control: the base calls as an alignment, and the chromatograms stacked (up to 48 clones); per-clone figures with --clone-figures
  alignedGroups(results).forEach((items, g) => { const tag = g ? `_${g + 1}` : ""; writeFigure(path.join(out, `clones_alignment${tag}.svg`), clonesAlignmentFigure(items)); if (items.length <= 48) writeFigure(path.join(out, `clones_aligned${tag}.svg`), clonesAlignedFigure(items)); });
  if (args["clone-figures"]) { const dir = path.join(out, "chromatograms"); mkdirSync(dir, { recursive: true }); results.forEach((r) => { if (r.ok && r.view) writeFigure(path.join(dir, `${String(r.name).replace(/[^\w.-]+/g, "_")}.svg`), cloneAlignedFigure(r.view, { name: r.name, decision: r.decision ? r.decision.decision : "", headline: headline(workflowKey, r) })); }); }
  for (const [file, svg] of pdfJobs) writeFileSync(file, await svgToPdf(svg));
  writeFileSync(path.join(out, "run.json"), JSON.stringify({ tool: "Assured QC", version: QC_VERSION, finished: new Date().toISOString(), seconds, threads, samples: results.length, samplesPerHour: Math.round((results.length / seconds) * 3600), arguments: { ...args, _: undefined }, geneRecord: genbank ? { file: path.basename(args.genbank), sha256: createHash("sha256").update(genbank).digest("hex"), options: codingOptions } : (designSpec && designSpec.coding && tasks.some((t) => t.options.codingInput) ? { file: "coding block of the design file", sha256: createHash("sha256").update(JSON.stringify(designSpec.coding)).digest("hex"), options: codingOptions } : null), ice: iceMeta || undefined, workflow: workflowKey, nuclease, pairingProblems: problems, decisions: tally(results.map((r) => r.decision)) }, null, 1));
  const t = tally(results.map((r) => r.decision)); console.log(`${results.length} samples in ${seconds.toFixed(1)} s on ${threads} thread(s) (${Math.round((results.length / seconds) * 3600)} per hour): accept ${t.accept}, hold ${t.hold}, re-sequence ${t["re-sequence"]}, reject ${t.reject}. Results in ${out}`);
  if (iceRows) { const s = iceSummary(iceRows); console.log(`ICE comparison (ice_comparison.csv): ${s.agree} agree, ${s.ice_fit_poor} with a poor ICE fit (not comparable), ${s.differ} differ${s.ice_failed ? `, ${s.ice_failed} failed` : ""}${s.not_run ? `, ${s.not_run} not run` : ""}.`); iceRows.filter((r) => r.status === "differ" || r.status === "ice_fit_poor").forEach((r) => console.log(`  ${r.sample}: Assured QC ${r.app}, ICE ${r.ice}${r.iceR2 !== null ? ` (R\u00b2 ${r.iceR2.toFixed(2)})` : ""}: ${r.note}`)); }
  const panel = panelOf(results); if (panel.length) console.log(`Take forward first: ${panel.join(", ")} (ranked in recommendation.csv).`);
  problems.forEach((p) => console.log(`  note: ${p}`));
  return 0;
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) main(process.argv.slice(2)).then((code) => process.exit(code));
