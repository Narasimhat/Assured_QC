// Runs Synthego ICE (a source checkout, not bundled) on Assured QC jobs, one Python process per share of the jobs.
// Separate processes, not multiprocessing: a process pool does not start on Windows. Nothing is uploaded; ICE runs on this machine.
import { spawn as nodeSpawn, spawnSync } from "node:child_process";
import { mkdirSync, writeFileSync, readFileSync, existsSync, rmSync } from "node:fs";
import { cpus } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { parseIceJson } from "../src/iceCompare.js";

const here = path.dirname(fileURLToPath(import.meta.url));
export const ICE_SCRIPT = path.join(here, "..", "scripts", "run_ice_jobs.py");

/** The Python to use: --ice-python, then $AQC_ICE_PYTHON, then python (python3 outside Windows). */
export const icePython = (args = {}) => args["ice-python"] || process.env.AQC_ICE_PYTHON || (process.platform === "win32" ? "python" : "python3");

/** ICE version from the checkout, and whether this Python has what ICE needs. Returns { ok, version, python, error }. */
export function checkIce(iceDir, python, run = spawnSync) {
  if (!iceDir || !existsSync(path.join(iceDir, "ice", "analysis.py"))) return { ok: false, error: `${iceDir || "(no folder)"} is not a Synthego ICE checkout (ice/analysis.py not found). Download https://github.com/synthego-open/ice and give its folder to --ice.` };
  let version = ""; try { version = (/__version__\s*=\s*['"]([^'"]+)['"]/.exec(readFileSync(path.join(iceDir, "ice", "__version__.py"), "utf8")) || [])[1] || ""; } catch { /* the version file is optional */ }
  const probe = run(python, ["-c", "import sys, numpy, scipy, pandas, Bio; print(sys.version.split()[0])"], { encoding: "utf8", timeout: 60000, stdio: ["ignore", "pipe", "pipe"] });
  if (probe.error) return { ok: false, version, error: `could not start ${python} (${probe.error.code || probe.error.message}). Use --ice-python with the path to a Python that has numpy, scipy, pandas and biopython.` };
  if (probe.status !== 0) return { ok: false, version, error: `${python} cannot import numpy, scipy, pandas and biopython, which ICE needs${probe.stderr ? ` (${String(probe.stderr).trim().split("\n").pop()})` : ""}. Use --ice-python with a Python that has them.` };
  return { ok: true, version, python: String(probe.stdout).trim() };
}

/**
 * Run the jobs. jobs: [{ control, edited, guides, donor }]. Returns a Map from job position to the script's result
 * ({ ok, summary, contribs } or { ok: false, error }). `spawnFn` is injectable for tests.
 */
export async function runIceJobs({ iceDir, python, jobs, workdir, concurrency, spawnFn = nodeSpawn, onProgress = () => {} }) {
  mkdirSync(workdir, { recursive: true });
  const n = Math.max(1, Math.min(concurrency || Math.max(1, cpus().length - 1), jobs.length || 1));
  const shards = Array.from({ length: n }, () => []); jobs.forEach((job, i) => shards[i % n].push({ i, job }));
  const results = new Map(); let finished = 0;
  await Promise.all(shards.map((shard, s) => new Promise((resolve) => {
    if (!shard.length) return resolve();
    const jobsFile = path.join(workdir, `jobs_${s}.json`), outFile = path.join(workdir, `out_${s}.json`);
    writeFileSync(jobsFile, JSON.stringify(shard.map(({ job }) => ({ control: job.control, edited: job.edited, guides: job.guides, donor: job.donor || "" }))));
    let stderr = ""; const fail = (error) => { shard.forEach(({ i }) => results.set(i, { ok: false, error })); finished += shard.length; onProgress(finished, jobs.length); resolve(); };
    let child; try { child = spawnFn(python, [ICE_SCRIPT, iceDir, jobsFile, outFile, "1"], { stdio: ["ignore", "pipe", "pipe"] }); } catch (e) { return fail(`could not start ${python}: ${e.message}`); }
    if (child.stdout && child.stdout.on) child.stdout.on("data", () => {});
    if (child.stderr && child.stderr.on) child.stderr.on("data", (d) => { stderr += String(d); });
    child.on("error", (e) => fail(`could not start ${python}: ${e.message}`));
    child.on("close", (code) => {
      let parsed = null; try { parsed = parseIceJson(readFileSync(outFile, "utf8")); } catch { /* handled below */ }
      if (!parsed) return fail(`ICE process ended with code ${code} and no result${stderr ? `: ${stderr.trim().split("\n").pop()}` : ""}`);
      shard.forEach(({ i }, k) => results.set(i, parsed[String(k)] || { ok: false, error: "ICE gave no result for this job" })); finished += shard.length; onProgress(finished, jobs.length); resolve();
    });
  })));
  try { rmSync(workdir, { recursive: true, force: true }); } catch { /* temporary files */ }
  return results;
}
