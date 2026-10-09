
import sys, os, json, time, traceback, tempfile, shutil, io, contextlib
ICE = sys.argv[1]; jobs_path = sys.argv[2]; out_path = sys.argv[3]; nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 1
sys.path.insert(0, ICE)

def run_one(args):
    i, job = args
    tmp = tempfile.mkdtemp(prefix="ice_", dir=os.path.dirname(out_path))
    t0 = time.time()
    try:
        from ice.analysis import single_sanger_analysis
        from io import StringIO
        from Bio import AlignIO
        from ice.classes import pair_alignment
        def to_clustal(self, aln, name1, name2):
            txt = ">" + name1 + "\n" + aln[0] + "\n>" + name2 + "\n" + aln[1]
            obj = list(AlignIO.parse(StringIO(txt), "fasta"))[0]
            b = StringIO(); AlignIO.write(obj, b, "clustal"); return b.getvalue().split("\n", 2)[2]
        pair_alignment.PairAlignment.align_list_to_clustal = to_clustal
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            res = single_sanger_analysis(job["control"], job["edited"], os.path.join(tmp, "res"), ",".join(job["guides"]), donor=(job["donor"] or None))
        if isinstance(res, str): res = json.loads(res)
        # keep the summary numbers and the top contributions only
        keep = {k: res.get(k) for k in ("ice", "ice_d", "rsq", "hdr_pct", "ko_score", "mean_discord_before", "mean_discord_after", "status", "notes", "guides", "warnings") if k in res}
        contribs = []
        cp = os.path.join(tmp, "res.contribs.json")
        if os.path.exists(cp):
            cj = json.load(open(cp))
            contribs = [dict(rel=c.get("rel_abundance"), seq=c.get("human_readable"), indel=(c.get("indel") or {}).get("total"), details=(c.get("indel") or {}).get("details"), wt=c.get("wt")) for c in (cj.get("contribs_list") or [])[:12]]
        return i, dict(ok=True, summary=keep, contribs=contribs, keys=list(res.keys()), seconds=round(time.time() - t0, 1))
    except Exception as e:
        return i, dict(ok=False, error=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-600:], seconds=round(time.time() - t0, 1))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
    jobs = json.load(open(jobs_path)); items = list(enumerate(jobs)); out = {}
    if nproc <= 1:
        for it in items:
            i, r = run_one(it); out[i] = r; print(i, r.get("ok"), r.get("seconds"), flush=True)
    else:
        import multiprocessing as mp
        with mp.Pool(nproc) as pool:
            for i, r in pool.imap_unordered(run_one, items):
                out[i] = r; print(i, r.get("ok"), r.get("seconds"), flush=True)
                json.dump(out, open(out_path, "w"))
    json.dump(out, open(out_path, "w"))
