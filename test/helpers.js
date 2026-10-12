import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";

const here = path.dirname(fileURLToPath(import.meta.url));
export const designDir = path.join(here, "..", "fixtures", "design");
export const simDir = path.join(here, "..", "fixtures", "sim");
export const loadSpec = (name) => JSON.parse(readFileSync(path.join(designDir, `${name}.design.json`), "utf8"));
export const loadTruth = () => JSON.parse(readFileSync(path.join(simDir, "truth.json"), "utf8"));

import { parseAbif } from "../src/abif.js";
import { reverseComplement } from "../src/seq.js";

const traceCache = new Map();
export function loadTrace(name) {
  if (!traceCache.has(name)) traceCache.set(name, parseAbif(readFileSync(path.join(simDir, name)), name));
  return traceCache.get(name);
}

/** The same trace as if it had been read from the opposite primer: samples, calls, quality and channels reversed and complemented. */
export function reverseTrace(trace) {
  const swap = { A: "T", C: "G", G: "C", T: "A" };
  const channels = {};
  for (const base of ["A", "C", "G", "T"]) channels[base] = Int16Array.from(trace.channels[swap[base]]).reverse();
  const last = trace.samples - 1;
  return {
    ...trace, name: `${trace.name}.rev`, channels,
    calls: reverseComplement(trace.calls),
    peakLocations: trace.peakLocations.map((p) => last - p).reverse(),
    quality: Uint8Array.from(trace.quality).reverse(),
  };
}
