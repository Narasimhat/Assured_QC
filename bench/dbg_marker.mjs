
import { readFileSync } from "node:fs";
import { runJobs } from "../src/runJobs.js";
const jobs = JSON.parse(readFileSync(process.argv[2],"utf8")); const idx = Number(process.argv[3]); const job = jobs[idx];
const file = (p, name) => { const b = readFileSync(p); return { name, buffer: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) }; };
const r = runJobs({ workflow: "snp", manual: { guides: job.guides.join(", "), donor: job.donor, gene: job.project }, sampleType: "clone", control: file(job.control,"c.ab1"), samples: [file(job.edited,"s.ab1")] });
const s = r.results[0]; console.log(Object.keys(s).join(","));
console.log(JSON.stringify(s.markers)); console.log(JSON.stringify(s.summary)); console.log(JSON.stringify(s.intervals));
