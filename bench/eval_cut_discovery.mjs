import { readFileSync, writeFileSync } from "node:fs";
import { parseAbif } from "../src/abif.js";
import { discoverCut } from "../src/cutDiscovery.js";
const jobs = JSON.parse(readFileSync(process.argv[2], "utf8")); const out = [];
const tr = (p) => { const b = readFileSync(p); return parseAbif(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength), p); };
jobs.forEach((j, i) => { try { const r = discoverCut({ control: tr(j.control), edited: tr(j.edited) }); out.push({ i, ...r }); } catch (e) { out.push({ i, ok: false, error: String(e).slice(0, 80) }); } });
writeFileSync(process.argv[3], JSON.stringify(out));
