"""Adapters that run one analysis tool on a list of jobs and return results in one common shape.

Job (JSON): {project, label, control, edited, guides: [spacer, ...], donor: "", sample_type: "clone"|"pool", kind: "KO"|"SNP", truth: {...} (optional)}
Result (per job, per tool): {ok, error, edited_pct, intended_pct, ko_pct, alleles: [{label, size, pct}], r2, seconds}
"""
import json, os, re, subprocess, sys, time, importlib.util

def _pct(x): return None if x is None else round(100 * float(x), 2)
def _size_from_label(label):
    m = re.search(r"([\-\u2212+])\s*(\d+)", str(label))
    if not m: return None
    return (-1 if m.group(1) in "-\u2212" else 1) * int(m.group(2))

def run_assured_qc(jobs, cfg, workdir):
    jin, jout = os.path.join(workdir, "aq_in.json"), os.path.join(workdir, "aq_out.json")
    json.dump(jobs, open(jin, "w"))
    root = cfg["assured_qc_dir"]
    p = subprocess.run([cfg["node"], "--preserve-symlinks", "--preserve-symlinks-main", os.path.join(root, "scripts", "run_validation.mjs"), jin, jout], cwd=root, capture_output=True, text=True, timeout=7200)
    if p.returncode: raise RuntimeError(p.stderr[-400:])
    out = []
    for r in json.load(open(jout)):
        if not r["ok"]: out.append(dict(ok=False, error=r.get("error"), seconds=r.get("seconds"))); continue
        s = r["summary"]
        out.append(dict(ok=True, error=None, edited_pct=s["editedPct"], intended_pct=s["intendedEditPct"] if jobs[r["index"]].get("donor") else None, ko_pct=s["koScorePct"], unexplained_pct=s.get("unexplainedPct"),
            alleles=[dict(label=c["label"], size=c.get("size"), pct=round(100 * c["fraction"], 2)) for c in r["contributions"]], r2=r["quality"]["r2"], seconds=r["seconds"],
            genotype=(r.get("genotype") or {}).get("category")))
    return out

def run_peaksplit(jobs, cfg, workdir):
    jin, jout = os.path.join(workdir, "ps_in.json"), os.path.join(workdir, "ps_out.json")
    json.dump(jobs, open(jin, "w"))
    p = subprocess.run([cfg["node"], os.path.join(os.path.dirname(__file__), "run_peaksplit.cjs"), cfg["peaksplit_dir"], jin, jout], capture_output=True, text=True, timeout=7200)
    if p.returncode: raise RuntimeError(p.stderr[-400:])
    out = []
    for r in json.load(open(jout)):
        if not r["ok"]: out.append(dict(ok=False, error=r.get("error"), seconds=r.get("seconds"))); continue
        out.append(dict(ok=True, error=None, edited_pct=_pct(1 - r["wt"]), intended_pct=_pct(r["ki"]), ko_pct=_pct(r["ko"]),
            alleles=[dict(label=t[0], size=_size_from_label(t[0]), pct=round(100 * t[2], 2)) for t in r["top"]], r2=r["r2"], seconds=r["seconds"], warnings=r.get("warnings")))
    return out

_te = None
def _load_traceedit(cfg):
    global _te
    if _te is None:
        spec = importlib.util.spec_from_file_location("te_engine", os.path.join(cfg["traceedit_dir"], "engine.py")); _te = importlib.util.module_from_spec(spec); spec.loader.exec_module(_te)
    return _te

def run_traceedit(jobs, cfg, workdir):
    te = _load_traceedit(cfg); out = []
    for j in jobs:
        t0 = time.time()
        try:
            c, s = te.read_ab1(j["control"], "control"), te.read_ab1(j["edited"], "sample")
            variants = te.donor_variants(c.sequence, j["donor"])["variants"] if j.get("donor") else None
            r = te.analyze(c, s, guides=j["guides"], variants=variants, max_deletion=40, max_insertion=2)
            outs = r["outcomes"]; nv = len(variants or [])
            wt = sum(o["fraction"] for o in outs if o["kind"] == "wild_type")
            full = sum(o["fraction"] for o in outs if o["kind"] == "substitution" and nv and o["label"].count("+") == nv - 1)
            ko = sum(o["fraction"] for o in outs if o["net_bp"] % 3 != 0)
            out.append(dict(ok=True, error=None, edited_pct=_pct(1 - wt), intended_pct=_pct(full) if nv else None, ko_pct=_pct(ko),
                alleles=[dict(label=o["label"], size=o["net_bp"], pct=round(100 * o["fraction"], 2)) for o in outs[:8]], r2=r["metrics"]["r_squared"], seconds=round(time.time() - t0, 2), warnings=r["warnings"]))
        except Exception as ex: out.append(dict(ok=False, error=str(ex)[:200], seconds=round(time.time() - t0, 2)))
    return out

def run_ice(jobs, cfg, workdir, nproc=8):
    slices = [(k, jobs[k::nproc]) for k in range(nproc)]; procs = []
    for k, sl in slices:
        if not sl: continue
        jin, jout = os.path.join(workdir, f"ice_in{k}.json"), os.path.join(workdir, f"ice_out{k}.json")
        json.dump(sl, open(jin, "w"))
        procs.append((k, sl, jout, subprocess.Popen([cfg["ice_python"], os.path.join(os.path.dirname(__file__), "..", "scripts", "run_ice_jobs.py"), cfg["ice_dir"], jin, jout, "1"],
                         stdout=open(os.path.join(workdir, f"ice_log{k}.txt"), "w"), stderr=subprocess.STDOUT, cwd=workdir)))
    results = [None] * len(jobs)
    for k, sl, jout, pr in procs:
        pr.wait(timeout=14400); part = json.load(open(jout))
        for n in range(len(sl)):
            v = part[str(n)]; idx = k + n * nproc
            if not v["ok"] or v["summary"].get("ice") is None: results[idx] = dict(ok=False, error=v.get("error") or v["summary"].get("notes") or "no result", seconds=v.get("seconds")); continue
            s = v["summary"]
            results[idx] = dict(ok=True, error=None, edited_pct=s["ice"], intended_pct=s["hdr_pct"] if jobs[idx].get("donor") else None, ko_pct=s["ko_score"],
                alleles=[dict(label=str(c["indel"]), size=c["indel"] if isinstance(c["indel"], int) else None, pct=round(100 * c["rel"], 2)) for c in v["contribs"][:8]], r2=s["rsq"], seconds=v.get("seconds"))
    return results

def run_tracy(jobs, cfg, workdir):
    raise NotImplementedError("Tracy adapter: see bench/README.md")

ADAPTERS = {"assured_qc": run_assured_qc, "peaksplit": run_peaksplit, "traceedit": run_traceedit, "ice": run_ice}
