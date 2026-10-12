import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { classifyFiles, buildPairs, readDesign, planTasks, runPool, inferWorkflow } from "../src/quick.js";
import { parseAbif } from "../src/abif.js";
import { analyseSample } from "../src/batch.js";
import { decide } from "../src/decision.js";

const here = path.dirname(fileURLToPath(import.meta.url)); const sim = (n) => path.join(here, "..", "fixtures", "sim", n);
const buffers = { "apoe_WT.700.ab1": "apoe_control.ab1", "apoe_A01_hom.700.ab1": "apoe_snp_clone_hom_edited.ab1", "apoe_A02_het.700.ab1": "apoe_snp_clone_het_edited.ab1", "apoe_B01_x.700.ab1": "nr2f2_control.ab1" };
const load = (name) => { const b = readFileSync(sim(buffers[name])); return { name, buffer: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) }; };

test("files are sorted into traces, sheets and designs", () => {
  const c = classifyFiles(["a.ab1", "B.ABI", "sheet.csv", "design.json", "report.html", "readme.pdf"]);
  assert.deepEqual([c.traces.length, c.sheets.length, c.designs.length, c.other.length], [2, 1, 2, 1]);
});

test("design: a JSON spec, an HTML report, and nothing", () => {
  const spec = JSON.parse(readFileSync(path.join(here, "..", "fixtures", "design", "apoe_r176c_snp.design.json"), "utf8"));
  assert.equal(readDesign(JSON.stringify(spec)).kind, "spec"); assert.equal(readDesign("").kind, "none");
  assert.throws(() => readDesign("<p>no</p>"), /No supported/);
  assert.equal(inferWorkflow(spec), "snp"); assert.equal(inferWorkflow(null), "knockout");
});

test("drop, pair, plan, run on a pool, decide: the whole flow without a browser", async () => {
  const names = ["apoe_WT.700.ab1", "apoe_A01_hom.700.ab1", "apoe_A02_het.700.ab1"];
  const { pairs, problems } = buildPairs({ traceNames: names }); assert.equal(pairs.length, 2); assert.equal(problems.length, 0); assert.ok(pairs.every((p) => p.confidence === "high"));
  const spec = JSON.parse(readFileSync(path.join(here, "..", "fixtures", "design", "apoe_r176c_snp.design.json"), "utf8"));
  const traces = new Map(names.map((n) => { const f = load(n); return [n, parseAbif(f.buffer, f.name)]; }));
  const { tasks, problems: planProblems } = planTasks({ pairs, getTrace: (n) => traces.get(n), design: { kind: "spec", spec } });
  assert.equal(planProblems.length, 0); assert.equal(tasks.length, 2); assert.ok(tasks.every((t) => t.workflow === "snp"));
  const progress = [];
  const spawn = () => ({ onmessage: null, terminate() {}, postMessage(message) { const task = message; setTimeout(() => {
    const controls = task.controls.map((f) => parseAbif(f.buffer, f.name)); const edited = task.edited.map((f) => parseAbif(f.buffer, f.name));
    const result = analyseSample({ name: task.name, spec: task.spec, control: controls[0], edited: edited[0], sampleType: task.sampleType });
    this.onmessage({ data: { index: task.index, result: { ...result, decision: decide(task.workflow, result, task.decisionOptions) } } });
  }, 0); } });
  const results = await runPool({ tasks, spawn, size: 2, payload: (task) => ({ ...task, controls: task.controls.map(load), edited: task.edited.map(load) }), onProgress: (p) => progress.push(p.done) });
  assert.deepEqual(progress.sort(), [1, 2]); const by = Object.fromEntries(results.map((r) => [r.name, r]));
  assert.equal(by["apoe_A01_hom.700"].decision.decision, "accept"); assert.equal(by["apoe_A02_het.700"].decision.decision, "hold");
});

test("a worker that dies is a failed sample, not a failed run", async () => {
  const tasks = [{ index: 0, name: "s1" }, { index: 1, name: "s2" }];
  const spawn = () => ({ onmessage: null, onerror: null, terminate() {}, postMessage() { setTimeout(() => this.onerror({ message: "boom" }), 0); } });
  const results = await runPool({ tasks, spawn, size: 1, payload: (t) => t });
  assert.equal(results.length, 2); assert.ok(results.every((r) => !r.ok && r.decision.decision === "re-sequence"));
});
