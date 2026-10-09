import { readFileSync, writeFileSync } from "node:fs";
import { autoPair } from "../src/sheet.js";
const data = JSON.parse(readFileSync(process.argv[2], "utf8")); const out = {};
for (const [proj, v] of Object.entries(data)) { const files = [...v.controls, ...Object.keys(v.pairs)]; const r = autoPair(files); const got = {}; r.pairs.forEach((p) => p.edited.forEach((e) => { got[e] = { control: p.control[0], confidence: p.confidence }; })); out[proj] = { got, unpaired: r.unpaired.map((u) => u.name), problems: r.problems }; }
writeFileSync(process.argv[3], JSON.stringify(out));
