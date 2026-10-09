// Worker thread: one sample per message. Reads its own files, so that nothing large crosses the thread boundary.
import { parentPort } from "node:worker_threads";
import { readFileSync } from "node:fs";
import { createHash } from "node:crypto";
import { parseAbif } from "../src/abif.js";
import { analyseSample } from "../src/batch.js";
import { decide } from "../src/decision.js";

const trace = (path) => { const bytes = readFileSync(path); return { trace: parseAbif(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), path.replace(/^.*[\\/]/, "")), sha256: createHash("sha256").update(bytes).digest("hex") }; };

parentPort.on("message", (task) => {
  const started = Date.now();
  try {
    const controls = task.controls.map(trace); const edited = task.edited.map(trace);
    const result = analyseSample({ name: task.name, spec: task.spec, control: controls[0].trace, edited: edited[0].trace, extraReads: edited.slice(1).map((e, i) => ({ name: `read ${i + 2}`, control: (controls[i + 1] || controls[0]).trace, edited: e.trace })), sampleType: task.sampleType, options: task.options || {} });
    const decision = decide(task.workflow, result, task.decisionOptions || {});
    parentPort.postMessage({ index: task.index, result: { ...result, decision, sampleType: task.sampleType, inputs: { controls: task.controls.map((p, i) => ({ file: p.replace(/^.*[\\/]/, ""), sha256: controls[i].sha256 })), edited: task.edited.map((p, i) => ({ file: p.replace(/^.*[\\/]/, ""), sha256: edited[i].sha256 })) }, seconds: (Date.now() - started) / 1000 } });
  } catch (error) {
    parentPort.postMessage({ index: task.index, result: { name: task.name, kind: "trace-pair", ok: false, error: String(error.message || error), warnings: [], decision: decide(task.workflow, { ok: false, error: String(error.message || error) }), seconds: (Date.now() - started) / 1000 } });
  }
});
