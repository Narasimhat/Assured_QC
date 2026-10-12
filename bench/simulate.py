"""Generate the labelled in-silico benchmark: forward and reverse reads, exact truth, one jobs.json.
Usage: python bench/simulate.py jobs_real.json outdir [--seed 20261009]
Controls and guides/donors come from a real jobs file (one entry per project is enough); every control is a real .ab1.
"""
import os, sys, json, argparse, hashlib, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import simlib

COMP = str.maketrans("ACGTN", "TGCAN")
def reverse_ops(ops, L):
    """The same edits in the coordinates of a read from the other end of an amplicon of L bases."""
    out = []
    for o in ops:
        if o[0] == "del": out.append(("del", L - o[1] - o[2], o[2]))
        elif o[0] == "ins": out.append(("ins", L - o[1], o[2][::-1].translate(COMP)))
        else: out.append(("sub", L - 1 - o[1], o[2].translate(COMP)))
    return out

def A(label, kind, size, pct, ops): return dict(label=label, kind=kind, size=size, pct=pct, ops=ops)
def dele(c, n): return [("del", c - n // 2, n)]
def feasible(ctl, ops_list):
    for ops in ops_list:
        for o in ops:
            if o[0] == "del" and (o[1] < 70 or ctl.usable_end - (o[1] + o[2]) < 100): return False
            if o[0] != "del" and (o[1] < 70 or ctl.usable_end - o[1] < 100): return False
    return True

def ko_scenarios(ctl, cuts, rng):
    a, b = sorted(cuts); span = b - a; S = {}
    between = [("del", a, span)] if span >= 10 else [("del", a - 10, span + 20)]
    bsize = -(between[0][2])
    W = lambda p: A("Wild type", "wt", 0, p, [])
    BT = lambda p, ops=None, size=None: A(f"Deletion between cuts ({-(size or bsize)} bp)", "between", size or bsize, p, ops or between)
    D = lambda n, c, p: A(f"-{n}", "indel", -n, p, dele(c, n))
    I = lambda seq, c, p: A(f"+{len(seq)}", "indel", len(seq), p, [("ins", c, seq)])
    base = lambda: "".join(rng.choice(list("ACGT"), 1))
    for r in (1, 2, 3): S[f"wt_rep{r}"] = [W(100)]
    for n, c in ((1, a), (2, b), (4, a), (7, b)): S[f"het_del{n}"] = [W(50), D(n, c, 50)]
    for n, c in ((1, a), (5, b)): S[f"hom_del{n}"] = [D(n, c, 100)]
    S["het_ins1"] = [W(50), I(base(), b, 50)]; S["hom_ins1"] = [I(base(), a, 100)]
    S["het_ins6"] = [W(50), I("".join(rng.choice(list("ACGT"), 6)), a, 50)]
    S["cmp_del3_ins1"] = [D(3, a, 50), I(base(), b, 50)]; S["cmp_del11_del5"] = [D(11, a, 50), D(5, b, 50)]
    for f in (2, 5, 10, 25, 50, 75, 100): S[f"between_f{f}"] = ([W(100 - f)] if f < 100 else []) + [BT(f)]
    j = [("del", a + 2, span - 3)] if span >= 14 else between
    S["between_jit_hom"] = [BT(100, j, -j[0][2])]; S["between_jit_het"] = [W(50), BT(50, j, -j[0][2])]
    for n in (40, 80, 150): S[f"large_del{n}_het"] = [W(50), D(n, a, 50)]
    S["large_del80_hom"] = [D(80, b, 100)]
    S["between_ins"] = [BT(50), I(base(), a, 50)]; S["between_del_wt"] = [BT(40), D(2, b, 30), W(30)]
    for k in (1, 2):
        w = rng.dirichlet(np.ones(6)) * 100; al = [W(w[0])]
        cands = [D(int(rng.integers(1, 12)), c, w[i + 1]) if rng.random() < .6 else I("".join(rng.choice(list("ACGT"), int(rng.integers(1, 5)))), c, w[i + 1]) for i, c in enumerate([a, b, a, b, a])]
        cands[0] = BT(w[1]); S[f"pool{k}"] = al + cands
    for f in (2, 5, 10): S[f"lod_ins1_f{f}"] = [W(100 - f), I(base(), a, f)]
    S["hard_between25"] = [W(75), BT(25)]; S["hard_het_del2"] = [W(50), D(2, a, 50)]
    return S

def snp_scenarios(ctl, cut, vars_, rng):
    S = {}; vs = sorted(vars_, key=lambda v: abs(v[0] - cut))
    W = lambda p: A("Wild type", "wt", 0, p, [])
    def H(p, subset=None, extra=None, label=None, kind=None):
        use = vs if subset is None else subset; ops = [("sub", int(pos), base) for pos, base in use] + (extra or [])
        full = subset is None
        return A(label or ("Intended edit (donor)" if full else "Partial conversion"), kind or ("hdr" if (full and not extra) else ("hdr_indel" if extra else "hdr_partial")), 0, p, ops)
    D = lambda n, p: A(f"-{n}", "indel", -n, p, dele(cut, n))
    I = lambda seq, p: A(f"+{len(seq)}", "indel", len(seq), p, [("ins", cut, seq)])
    base = lambda: "".join(rng.choice(list("ACGT"), 1))
    for r in (1, 2, 3): S[f"wt_rep{r}"] = [W(100)]
    for f in (2, 5, 10, 25, 50, 75, 100): S[f"hdr_f{f}"] = ([W(100 - f)] if f < 100 else []) + [H(f)]
    S["hdr_het_ins"] = [H(50), I(base(), 50)]; S["hdr_het_del"] = [H(50), D(3, 50)]
    ins_other = [("ins", cut, base())]
    if not any(abs(p - cut) <= 1 for p, _ in vs): S["hdrindel_wt"] = [W(50), H(50, extra=ins_other, label="Donor allele with +1 indel", kind="hdr_indel")]; S["hdrindel_wt"][1]["size"] = 1
    if len(vs) >= 2:
        near = vs[:1]; far = vs[1:]
        S["partial_near_hom"] = [H(100, near)]; S["partial_far_het"] = [W(50), H(50, far)]
        S["partial_mix"] = [W(20), H(60, near), H(20)]
    S["het_del2"] = [W(50), D(2, 50)]; S["hom_ins1"] = [I(base(), 100)]; S["het_del8"] = [W(50), D(8, 50)]
    S["mix1"] = [H(40), I(base(), 30), W(30)]; S["mix2"] = [H(30), D(4, 30), W(20), I(base(), 20)]
    S["stress_hdr10"] = [W(90), H(10)]; S["stress_hdr50"] = [W(50), H(50)]
    S["hard_hdr25"] = [W(75), H(25)]; S["hard_het_del2"] = [W(50), D(2, 50)]
    return S

def main(jobs_real, outdir, seed=20261009):
    real = json.load(open(jobs_real)); os.makedirs(os.path.join(outdir, "traces"), exist_ok=True)
    seen, jobs = {}, []
    for j in real:
        seen.setdefault((j["project"], j["control"]), j)
    ctls = {}
    for n, ((proj, cpath), j) in enumerate(seen.items()):
        ctl = ctls[cpath] = simlib.Control(cpath); rng = np.random.default_rng(seed + n)
        guides, donor = j["guides"], j.get("donor") or ""; cuts = [ctl.find_guide(g) for g in guides]
        S = snp_scenarios(ctl, cuts[0], ctl.donor_variants(donor.upper()), rng) if donor else ko_scenarios(ctl, cuts, rng)
        ctag = os.path.splitext(os.path.basename(cpath))[0]
        rctl_path = os.path.join(outdir, "traces", f"{proj}_{ctag}_ctl_rev.ab1")
        amplicon = ctl.truncated(min(ctl.L, ctl.usable_end + 40)); rctl = amplicon.reversed()
        simlib.write_ab1(rctl_path, rctl.ch, rctl.calls, rctl.ploc, rctl.pcon, "ctl_rev")
        for name, alleles in S.items():
            ops_list = [a["ops"] for a in alleles]
            if not feasible(ctl, ops_list): continue
            tot = sum(a["pct"] for a in alleles); assert abs(tot - 100) < 1e-6, (name, tot)
            h = int(hashlib.sha1(f"{seed}|{proj}|{ctag}|{name}".encode()).hexdigest()[:8], 16); trng = np.random.default_rng(h)
            stress = name.startswith(("stress", "hard"))
            noise = 0.08 if name.startswith("hard") else (0.04 if stress else float(trng.choice([0.0, 0.01, 0.02])))
            dye = tuple(1 + trng.uniform(-0.1, 0.1, 4)) if (stress or trng.random() < 0.3) else (1, 1, 1, 1)
            amp = float(trng.uniform(0.7, 1.3))
            items = [(simlib.build_allele(ctl.L, a["ops"]), a["pct"] / 100) for a in alleles]
            ch, calls, ploc, qual = simlib.render(ctl, items, trng, noise=noise, dye=dye, amp=amp, jitter=0.08)
            stem = f"{proj}_{ctag}_{name}"; fpath = os.path.join(outdir, "traces", stem + "_F.ab1"); simlib.write_ab1(fpath, ch, calls, ploc, qual, stem)
            rpath = os.path.join(outdir, "traces", stem + "_R.ab1")
            try:
                ritems = [(simlib.build_allele(rctl.L, reverse_ops(a["ops"], amplicon.L)), a["pct"] / 100) for a in alleles]
                rrng = np.random.default_rng(h + 1)
                rch, rcalls, rploc, rq = simlib.render(rctl, ritems, rrng, noise=noise, dye=tuple(1 + rrng.uniform(-0.1, 0.1, 4)) if dye != (1, 1, 1, 1) else (1, 1, 1, 1), amp=float(rrng.uniform(0.7, 1.3)), jitter=0.08)
                simlib.write_ab1(rpath, rch, rcalls, rploc, rq, stem + "_R")
            except Exception as ex: rpath = None
            wt = sum(a["pct"] for a in alleles if a["kind"] == "wt")
            frame = lambda a: a["kind"] in ("indel", "between", "hdr_indel") and (a["size"] % 3 != 0 or abs(a["size"]) >= 21)
            truth = dict(edited_pct=round(100 - wt, 3), intended_pct=round(sum(a["pct"] for a in alleles if a["kind"] == "hdr"), 3) if donor else None,
                         ko_pct=round(sum(a["pct"] for a in alleles if frame(a)), 3), alleles=[{k: a[k] for k in ("label", "kind", "size", "pct")} for a in alleles])
            jobs.append(dict(project=proj, label=stem, control=cpath, edited=fpath, control_rev=rctl_path, edited_rev=rpath, guides=guides, donor=donor,
                             sample_type="pool" if name.startswith("pool") else "clone", kind="SNP" if donor else "KO", truth=truth,
                             sim=dict(scenario=name, seed=h, noise=noise, dye=[round(float(x), 3) for x in dye], amp=round(amp, 3), cuts=[int(c) for c in cuts])))
    json.dump(jobs, open(os.path.join(outdir, "jobs.json"), "w"))
    return jobs

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("jobs_real"); ap.add_argument("outdir"); ap.add_argument("--seed", type=int, default=20261009); a = ap.parse_args()
    jobs = main(a.jobs_real, a.outdir, a.seed); print(len(jobs), "labelled traces written to", a.outdir)
