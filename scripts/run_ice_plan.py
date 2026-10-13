"""Run the ICE jobs of an Assured QC run (ice_jobs.json from `--ice-plan`) and write ice_raw.json for `--ice-raw`.

    python scripts/run_ice_plan.py <ICE checkout> ice_jobs.json ice_raw.json [processes]

Works wherever Python and ICE do (Windows, Linux, a cluster node): the jobs are split over separate processes, because a process pool does not start on
Windows. The traces are read from the paths in ice_jobs.json and recorded by SHA-256, so Assured QC accepts the results only for the same traces, guides
and donors. Synthego ICE is used under its research and non-profit licence and is not redistributed.
"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def sha256(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def main(ice_dir, jobs_path, out_path, nproc=4):
    plan = json.load(open(jobs_path))
    jobs = plan["jobs"]
    work = tempfile.mkdtemp(prefix="ice_plan_", dir=os.path.dirname(os.path.abspath(out_path)))
    shards = [[] for _ in range(max(1, min(nproc, len(jobs) or 1)))]
    for i, job in enumerate(jobs):
        shards[i % len(shards)].append((i, job))
    running = []
    for s, shard in enumerate(shards):
        if not shard:
            continue
        jobs_file, out_file = os.path.join(work, "jobs_%d.json" % s), os.path.join(work, "out_%d.json" % s)
        json.dump([{"control": j["control"], "edited": j["edited"], "guides": j["guides"], "donor": j.get("donor") or ""} for _, j in shard], open(jobs_file, "w"))
        log_file = os.path.join(work, "log_%d.txt" % s)
        log = open(log_file, "w")
        running.append((subprocess.Popen([sys.executable, os.path.join(HERE, "run_ice_jobs.py"), ice_dir, jobs_file, out_file, "1"], stdout=log, stderr=subprocess.STDOUT), out_file, shard, log, log_file))
    results = {}
    for proc, out_file, shard, log, log_file in running:
        proc.wait()
        log.close()
        try:
            data = json.load(open(out_file))
        except Exception:
            data = {}
        tail = ""
        try:
            tail = open(log_file).read().strip().splitlines()[-1]
        except Exception:
            pass
        for k, (i, _) in enumerate(shard):
            results[i] = data.get(str(k)) or {"ok": False, "error": "ICE process ended without a result" + (": " + tail if tail else "")}
    shutil.rmtree(work, ignore_errors=True)
    version = ""
    try:
        version = re.search(r"__version__\s*=\s*['\"]([^'\"]+)", open(os.path.join(ice_dir, "ice", "__version__.py")).read()).group(1)
    except Exception:
        pass
    records = [{"sample": j["sample"], "donorName": j.get("donorName") or "", "guides": j["guides"], "donor": j.get("donor") or "", "control": os.path.basename(j["control"]), "controlSha256": sha256(j["control"]),
                "edited": os.path.basename(j["edited"]), "editedSha256": sha256(j["edited"]), "result": results[i]} for i, j in enumerate(jobs)]
    json.dump({"tool": "Synthego ICE", "version": version, "python": sys.version.split()[0], "records": records}, open(out_path, "w"), indent=1)
    print("%d ICE run(s), %d ok -> %s" % (len(records), sum(1 for r in records if r["result"].get("ok")), out_path))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 4)
