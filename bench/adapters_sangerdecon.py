"""sangerdecon adapter for the Assured_QC benchmark harness (kept in its own file so the frozen bench/adapters.py is untouched).

Result shape (bench/README.md): {ok, error, edited_pct, intended_pct, ko_pct, alleles: [{label, size, pct}], r2, seconds}

Two variants are written because the harness scores the point estimate against a 2-point detection threshold:
  sangerdecon           raw point estimate (1 - wild type), comparable to the other tools' headline number
  sangerdecon_reported  the same fit with sangerdecon's own detection limits applied (indels < 3% and donor changes < 10% reported as 0)

Usage: python adapters_sangerdecon.py jobs.json outdir [--sangerdecon-dir DIR] [--procs N]
"""
import argparse, json, os, sys, time


def _label(kind, net, top):
    if kind == "WT":
        return "Wild type"
    if kind == "SUB":
        return "Intended edit (donor)" if "(full)" in str(top) else "Partial conversion"
    if kind == "INV":
        return "Inversion"
    return f"{net:+d}"


def _run_one(args):
    job, sd_dir = args
    sys.path.insert(0, sd_dir)
    from sangerdecon import analyze as AN
    t0 = time.time()
    try:
        r = AN.analyze_pair(job["control"], job["edited"], list(job["guides"]), donor=(job.get("donor") or None))
    except Exception as ex:   # keep the batch going; the scorer counts it as not analysed
        return dict(ok=False, error=f"{type(ex).__name__}: {str(ex)[:150]}", seconds=round(time.time() - t0, 2))
    if not r.get("ok"):
        return dict(ok=False, error=",".join(r.get("flags", [])) or "not analysed", seconds=round(time.time() - t0, 2))
    has_donor = bool(job.get("donor"))
    alleles = [dict(label=_label(c.kind, int(c.net), c.top_allele), size=(int(c.net) if c.kind not in ("SUB", "WT") else 0), pct=round(100 * float(c.fraction), 2))
               for c in r["classes"].head(8).itertuples()]
    return dict(ok=True, error=None, edited_pct=round(100 * r["editing"], 3), edited_pct_reported=round(100 * r["editing_reported"], 3),
                intended_pct=(round(100 * r["donor_full"], 3) if has_donor else None), ko_pct=round(100 * r["ko_score"], 3),
                alleles=alleles, r2=float(r["r2"]), seconds=round(time.time() - t0, 2), flags=r.get("flags", []), genotype=r.get("genotype"))


def run_sangerdecon(jobs, sd_dir, procs=8, workdir="."):
    """Slices run in separate interpreter processes (a process pool is not available in every sandbox)."""
    import subprocess
    os.makedirs(workdir, exist_ok=True); running = []
    for k in range(procs):
        sl = list(range(k, len(jobs), procs))
        if not sl:
            continue
        jin, jout = os.path.join(workdir, f"sd_in{k}.json"), os.path.join(workdir, f"sd_out{k}.json")
        json.dump([jobs[i] for i in sl], open(jin, "w"))
        running.append((sl, jout, subprocess.Popen([sys.executable, os.path.abspath(__file__), "--worker", jin, jout, "--sangerdecon-dir", sd_dir],
                                                   stdout=open(os.path.join(workdir, f"sd_log{k}.txt"), "w"), stderr=subprocess.STDOUT)))
    out = [None] * len(jobs)
    for sl, jout, pr in running:
        pr.wait()
        res = json.load(open(jout)) if os.path.exists(jout) else [dict(ok=False, error="worker failed")] * len(sl)
        for i, r in zip(sl, res):
            out[i] = r
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs"); ap.add_argument("outdir")
    ap.add_argument("--sangerdecon-dir", default=os.environ.get("SANGERDECON_DIR", "."))
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--worker", action="store_true", help="internal: analyse a slice (jobs=in.json, outdir=out.json)")
    a = ap.parse_args()
    if a.worker:
        jobs = json.load(open(a.jobs)); sd = os.path.abspath(a.sangerdecon_dir)
        json.dump([_run_one((j, sd)) for j in jobs], open(a.outdir, "w")); sys.exit(0)
    jobs = json.load(open(a.jobs)); os.makedirs(a.outdir, exist_ok=True)
    t0 = time.time(); res = run_sangerdecon(jobs, os.path.abspath(a.sangerdecon_dir), a.procs, workdir=os.path.join(a.outdir, "sd_work"))
    json.dump(res, open(os.path.join(a.outdir, "sangerdecon.json"), "w"))
    rep = [dict(r, edited_pct=r.get("edited_pct_reported")) if r["ok"] else r for r in res]
    json.dump(rep, open(os.path.join(a.outdir, "sangerdecon_reported.json"), "w"))
    print(f"sangerdecon: {sum(r['ok'] for r in res)}/{len(res)} analysed in {time.time() - t0:.0f}s")
