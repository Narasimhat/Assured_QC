"""Scorecard: compare tools against truth on a labelled job file (simulated, known-mixture, or any job file with truth).
Usage: python bench/score.py jobs.json results_dir out_prefix [--split dev|test|all]
Definitions are fixed in docs/scorecard.md; do not change them after results are seen (hash in bench/FROZEN.json).
"""
import argparse, json, os, re, sys
import numpy as np, pandas as pd

DEV_PROJECTS = {"L01_KO_SNP_editing", "L02_KO_SNP_editing", "L04_KO_SNP_editing", "L07_SNP_editing", "L09_Exn20_SNP_editing", "L11_SNP_editing"}
DETECT_T = 2.0          # points of edited signal; a trace is "detected" at or above this
BANDS = [(0, 0, "0"), (0.0001, 5, "(0,5]"), (5.0001, 25, "(5,25]"), (25.0001, 75, "(25,75]"), (75.0001, 100, "(75,100]")]
FAMILIES = [("wt", r"^wt_rep"), ("hom_or_het_small_indel", r"^(het|hom)_(del|ins)"), ("compound_het", r"^cmp_"), ("two_cut_graded", r"^between_f"), ("two_cut_junction_shift", r"^between_jit"),
            ("two_cut_plus_indel", r"^between_(ins|del)"), ("large_deletion", r"^large_del"), ("pool", r"^pool"), ("low_fraction_indel", r"^lod_"), ("hdr_graded", r"^hdr_f"),
            ("hdr_with_indel", r"^(hdr_het|hdrindel)"), ("partial_conversion", r"^partial_"), ("indel_only_snp_locus", r"^(het_del2|hom_ins1|het_del8)$"), ("mixed_alleles", r"^mix"), ("noisy", r"^(stress|hard)")]
def family(name):
    for fam, rx in FAMILIES:
        if re.search(rx, name): return fam
    return "other"
def band(v):
    for lo, hi, lab in BANDS:
        if lo <= v <= hi: return lab
    return "?"

def tool_alleles_by_size(alleles):
    out = {}
    for a in alleles or []:
        if a.get("size") is None: continue
        out[a["size"]] = out.get(a["size"], 0) + a["pct"]
    return out

def load(jobs_path, rdir, split):
    jobs = json.load(open(jobs_path)); tools = [f[:-5] for f in sorted(os.listdir(rdir)) if f.endswith(".json") and not f.startswith(("aq_", "ps_", "ice_"))]
    keep = [i for i, j in enumerate(jobs) if split == "all" or (split == "dev") == (j["project"] in DEV_PROJECTS)]
    res = {t: json.load(open(os.path.join(rdir, t + ".json"))) for t in tools}
    return jobs, keep, res

def per_trace(jobs, keep, res):
    rows = []
    for t, R in res.items():
        for i in keep:
            j, r, tr = jobs[i], R[i], jobs[i]["truth"]
            ok = bool(r.get("ok")); sim = j.get("sim", {}); scn = sim.get("scenario", j.get("label"))
            est = (r.get("edited_pct") if ok else None); estI = r.get("intended_pct") if ok else None; estK = r.get("ko_pct") if ok else None
            rows.append(dict(tool=t, i=i, project=j["project"], scenario=scn, family=family(scn), kind=j["kind"], ok=ok, truth_edited=tr["edited_pct"], est_edited=est,
                err=None if est is None else est - tr["edited_pct"], truth_intended=tr.get("intended_pct"), est_intended=estI, truth_ko=tr.get("ko_pct"), est_ko=estK, seconds=r.get("seconds"),
                alleles=r.get("alleles") if ok else None, truth_alleles=tr["alleles"]))
    return pd.DataFrame(rows)

def allele_recovery(df):
    """A truth allele (indel or two-cut deletion, summed over alleles of equal net size) at >= 10 % is recovered if the tool reports that net size at between half and twice the truth fraction."""
    out = []
    for _, r in df.iterrows():
        truth = {}
        for a in r.truth_alleles:
            if a["kind"] in ("indel", "between"): truth[a["size"]] = truth.get(a["size"], 0) + a["pct"]
        mine = tool_alleles_by_size(r.alleles) if r.ok else {}
        for size, pct in truth.items():
            if pct < 10: continue
            got = mine.get(size, 0)
            out.append(dict(tool=r.tool, family=r.family, i=r.i, truth_pct=pct, est_pct=got, recovered=bool(r.ok and pct / 2 <= got <= pct * 2)))
    return pd.DataFrame(out)

def boot_ci(values, n=2000, seed=1):
    rng = np.random.default_rng(seed); v = np.asarray(values, float)
    if len(v) == 0: return (np.nan, np.nan)
    m = rng.choice(v, (n, len(v))).mean(axis=1); return tuple(np.percentile(m, [2.5, 97.5]))

def scorecard(df):
    rows = []
    for t, g in df.groupby("tool"):
        a = g[g.ok]; imp = g.assign(e=np.where(g.ok, g.est_edited, 0.0)); imp["abs_imputed"] = (imp.e - imp.truth_edited).abs()
        wt = g[g.family == "wt"]; fp = np.mean([(r.ok and r.est_edited >= DETECT_T) for r in wt.itertuples()]) if len(wt) else np.nan
        row = dict(tool=t, n=len(g), coverage=a.shape[0] / len(g), MAE_edited_analysed=a.err.abs().mean(), MAE_edited_all_imputed=imp.abs_imputed.mean(),
                   MAE_edited_ci_lo=boot_ci(a.err.abs())[0], MAE_edited_ci_hi=boot_ci(a.err.abs())[1], median_abs_err=a.err.abs().median(), p95_abs_err=a.err.abs().quantile(.95), FP_rate_WT=fp, seconds_median=g.seconds.median())
        for lab in [b[2] for b in BANDS]:
            x = a[[band(v) == lab for v in a.truth_edited]]; row[f"MAE_band_{lab}"] = x.err.abs().mean() if len(x) else np.nan
        sn = a[a.truth_intended.notna() & a.est_intended.notna()]; row["MAE_intended"] = (sn.est_intended - sn.truth_intended).abs().mean() if len(sn) else np.nan
        pc = a[(a.truth_intended.notna()) & (a.truth_intended < 5) & (a.family.isin(["partial_conversion"])) & a.est_intended.notna()]
        row["false_intended_rate_partial"] = (pc.est_intended > 10).mean() if len(pc) else np.nan
        kk = a[a.truth_ko.notna() & a.est_ko.notna()]; row["MAE_ko"] = (kk.est_ko - kk.truth_ko).abs().mean() if len(kk) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index("tool")

def lod_table(df):
    fam = df[df.family.isin(["two_cut_graded", "hdr_graded", "low_fraction_indel"])].copy()
    fam["level"] = fam.truth_edited.round(0); out = []
    for (t, f, lv), g in fam.groupby(["tool", "family", "level"]):
        out.append(dict(tool=t, family=f, level=lv, n=len(g), sensitivity=np.mean([(r.ok and r.est_edited >= DETECT_T) for r in g.itertuples()]), median_est=g.est_edited.median()))
    return pd.DataFrame(out)

def by_family(df):
    a = df[df.ok]
    return a.assign(abs_err=a.err.abs()).pivot_table(index="family", columns="tool", values="abs_err", aggfunc="mean").round(2).join(df.groupby(["family"]).size().rename("n_traces") // df.tool.nunique())

def paired_difference(df, ref="assured_qc", seed=2, n=2000):
    """Per-trace |err| difference (comparator - reference) over traces that both analysed; bootstrap 95% CI of the mean. Positive = reference is better."""
    rng = np.random.default_rng(seed); p = df.assign(a=df.err.abs()).pivot_table(index="i", columns="tool", values="a")
    if ref not in p.columns: return pd.DataFrame()
    out = []
    for t in p.columns:
        if t == ref: continue
        d = (p[t] - p[ref]).dropna().values
        m = rng.choice(d, (n, len(d))).mean(axis=1); out.append(dict(comparator=t, n_common=len(d), mean_diff=d.mean(), ci_lo=np.percentile(m, 2.5), ci_hi=np.percentile(m, 97.5)))
    return pd.DataFrame(out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("jobs"); ap.add_argument("results"); ap.add_argument("prefix"); ap.add_argument("--split", default="all"); ap.add_argument("--ref", default="assured_qc"); a = ap.parse_args()
    jobs, keep, res = load(a.jobs, a.results, a.split); df = per_trace(jobs, keep, res)
    sc = scorecard(df); ar = allele_recovery(df)
    sc["allele_recovery"] = ar.groupby("tool").recovered.mean(); sc["allele_recovery_large_del"] = ar[ar.family == "large_deletion"].groupby("tool").recovered.mean()
    sc.round(3).to_csv(a.prefix + "_scorecard.csv"); lod_table(df).round(3).to_csv(a.prefix + "_detection.csv", index=False); by_family(df).to_csv(a.prefix + "_by_family.csv"); paired_difference(df, a.ref).round(3).to_csv(a.prefix + "_paired.csv", index=False)
    df.drop(columns=["alleles", "truth_alleles"]).to_csv(a.prefix + "_per_trace.csv", index=False)
    print(sc.round(2).T.to_string())
