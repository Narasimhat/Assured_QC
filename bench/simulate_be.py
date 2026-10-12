"""Simulated base-editing traces (A>G or C>T in protospacer positions 4-8) from real controls, with exact per-position conversion truth.
Each position converts independently with its own probability, so a trace is a mixture of up to 2^k alleles.
Usage: python bench/simulate_be.py jobs_real.json outdir [--seed 20261011]"""
import os, sys, json, argparse, hashlib, itertools, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import simlib
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

def place(ctl, guide):
    for strand, s in (("+", guide), ("-", simlib.rc(guide))):
        i = ctl.calls.find(s)
        if i >= 0: return strand, i, i + len(s)
    raise ValueError("guide not found")

def main(jobs_real, outdir, seed):
    real = json.load(open(jobs_real)); os.makedirs(os.path.join(outdir, "traces"), exist_ok=True); seen = {}
    for j in real: seen.setdefault((j["project"], j["control"]), j)
    jobs = []
    for n, ((proj, cpath), j) in enumerate(seen.items()):
        ctl = simlib.Control(cpath); guide = j["guides"][0]; strand, start, end = place(ctl, guide)
        for editor, frm, to in (("ABE", "A", "G"), ("CBE", "C", "T")):
            pos = []
            for p in range(4, 9):
                idx = start + p - 1 if strand == "+" else end - p
                base = ctl.calls[idx]; gb = base if strand == "+" else COMP[base]
                if gb == frm: pos.append((p, idx, base, to if strand == "+" else COMP[to]))
            if not pos: continue
            for rep in range(3):
                h = int(hashlib.sha1(f"{seed}|{proj}|{cpath}|{editor}|{rep}".encode()).hexdigest()[:8], 16); rng = np.random.default_rng(h)
                rates = {p: float(np.round(rng.choice([0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95]), 2)) for p, *_ in pos}
                target = pos[int(rng.integers(len(pos)))][0]; rates[target] = float(np.round(rng.uniform(0.15, 0.95), 2))
                alleles = []
                for combo in itertools.product([0, 1], repeat=len(pos)):
                    prob = 1.0
                    for c, (p, *_ ) in zip(combo, pos): prob *= rates[p] if c else 1 - rates[p]
                    if prob > 1e-6: alleles.append((simlib.build_allele(ctl.L, [("sub", idx, alt) for c, (p, idx, ref, alt) in zip(combo, pos) if c]), prob))
                tot = sum(a for _, a in alleles); alleles = [(it, f / tot) for it, f in alleles]
                noise = float(rng.choice([0.0, 0.01, 0.02]))
                ch, calls, ploc, qual = simlib.render(ctl, alleles, rng, noise=noise, amp=float(rng.uniform(0.7, 1.3)), jitter=0.08)
                stem = f"{proj}_{os.path.splitext(os.path.basename(cpath))[0]}_{editor}_{rep}"; path = os.path.join(outdir, "traces", stem + ".ab1"); simlib.write_ab1(path, ch, calls, ploc, qual, stem)
                jobs.append(dict(project=proj, label=stem, control=cpath, edited=path, guides=[guide], editor=editor, target=target, sample_type="pool",
                                 truth=dict(conversion={str(p): rates[p] * 100 for p, *_ in pos}, positions={str(p): idx for p, idx, *_ in pos}), sim=dict(noise=noise, seed=h)))
    json.dump(jobs, open(os.path.join(outdir, "jobs.json"), "w")); return jobs

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("jobs_real"); ap.add_argument("outdir"); ap.add_argument("--seed", type=int, default=20261011); a = ap.parse_args()
    print(len(main(a.jobs_real, a.outdir, a.seed)), "base-editing traces")
