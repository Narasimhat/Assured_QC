"""Secondary tables for the scorecard that need the full Assured QC output (intervals, tiers, decisions).
Usage: python bench/extras.py jobs.json results_dir out_prefix [--split dev|test|all]
  results_dir contains aq_out.json (raw engine output, written by the benchmark adapter) and dec_out.json (bench/run_decisions.mjs).
Definitions (fixed before the held-out run, see docs/scorecard.md):
  interval coverage   share of truths inside the reported interval (edited share, intended-edit share when a donor is used, knockout score)
  detected flag       interval lower bound of the edited share >= 1 point; sensitivity at each graded level and false-positive rate on wild-type traces
  confidence tiers    mean and 95th-percentile absolute error of the edited share per tier
  decision            clones only; 'worthy' = truth intended-edit share (donor designs) or knockout share (others) >= 90 %;
                      false accept = accept on a clone that is not worthy; true accept = worthy clone accepted
"""
import argparse, json, os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from score import DEV_PROJECTS, family


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("jobs"); ap.add_argument("results"); ap.add_argument("prefix"); ap.add_argument("--split", default="test"); a = ap.parse_args()
    jobs = json.load(open(a.jobs)); raw = {r["index"]: r for r in json.load(open(os.path.join(a.results, "aq_out.json")))}
    dec = {d["index"]: d for d in json.load(open(os.path.join(a.results, "dec_out.json")))}
    keep = [i for i, j in enumerate(jobs) if a.split == "all" or (a.split == "dev") == (j["project"] in DEV_PROJECTS)]
    rows = []
    for i in keep:
        j = jobs[i]; r = raw.get(i); d = dec.get(i); tr = j["truth"]; scn = j["sim"]["scenario"]
        row = dict(i=i, project=j["project"], scenario=scn, family=family(scn), donor=bool(j.get("donor")), sample_type=j.get("sample_type", "clone"), truth_edited=tr["edited_pct"], truth_intended=tr.get("intended_pct"), truth_ko=tr.get("ko_pct"),
                   ok=bool(r and r["ok"]), decision=(d or {}).get("decision", {}).get("decision") if d else None, tier=(d or {}).get("tier"))
        if r and r["ok"]:
            iv = r["intervals"]; s = r["summary"]
            row.update(est_edited=s["editedPct"], e_lo=iv["editedPct"][0], e_hi=iv["editedPct"][1], i_lo=iv["intendedEditPct"][0], i_hi=iv["intendedEditPct"][1], k_lo=iv["koScorePct"][0], k_hi=iv["koScorePct"][1], detected=bool(iv["detected"]))
        rows.append(row)
    df = pd.DataFrame(rows); ok = df[df.ok].copy(); tol = 0.05
    out = {}
    cov = lambda lo, hi, t: float(((lo - tol <= t) & (t <= hi + tol)).mean())
    out["n_analysed"] = int(len(ok)); out["n_total"] = int(len(df))
    out["coverage_edited"] = cov(ok.e_lo, ok.e_hi, ok.truth_edited); out["width_edited_median"] = float((ok.e_hi - ok.e_lo).median())
    oi = ok[ok.donor & ok.truth_intended.notna()]; out["coverage_intended"] = cov(oi.i_lo, oi.i_hi, oi.truth_intended); out["width_intended_median"] = float((oi.i_hi - oi.i_lo).median()); out["n_intended"] = int(len(oi))
    ok_k = ok[ok.truth_ko.notna()]; out["coverage_ko"] = cov(ok_k.k_lo, ok_k.k_hi, ok_k.truth_ko); out["width_ko_median"] = float((ok_k.k_hi - ok_k.k_lo).median())
    # coverage by family (edited share)
    byfam = ok.groupby("family").apply(lambda g: pd.Series(dict(n=len(g), coverage_edited=cov(g.e_lo, g.e_hi, g.truth_edited), width_median=float((g.e_hi - g.e_lo).median()))), include_groups=False).round(3)
    byfam.to_csv(a.prefix + "_interval_by_family.csv")
    # detected flag
    wt = ok[ok.family == "wt"]; out["detected_flag_FP_wt"] = float(wt.detected.mean()) if len(wt) else None; out["n_wt"] = int(len(wt))
    lod = ok[ok.family.isin(["two_cut_graded", "hdr_graded", "low_fraction_indel"])].copy(); lod["level"] = lod.truth_edited.round(0)
    lod.groupby(["family", "level"]).apply(lambda g: pd.Series(dict(n=len(g), sensitivity_detected_flag=float(g.detected.mean()), sensitivity_point_ge2=float((g.est_edited >= 2).mean()))), include_groups=False).round(3).to_csv(a.prefix + "_lod_detected_flag.csv")
    # tiers
    ok["abs_err"] = (ok.est_edited - ok.truth_edited).abs()
    tiers = ok.dropna(subset=["tier"]).groupby("tier").abs_err.agg(n="size", mean="mean", p95=lambda x: x.quantile(.95), share_gt5=lambda x: (x > 5).mean()).round(3); tiers["share_of_traces"] = (tiers.n / tiers.n.sum()).round(3); tiers.to_csv(a.prefix + "_tiers.csv")
    # decisions (clones)
    cl = df[(df.sample_type == "clone") & df.decision.notna()].copy(); cl["truth_goal"] = np.where(cl.donor, cl.truth_intended, cl.truth_ko); cl["worthy"] = cl.truth_goal >= 90
    out["clones"] = int(len(cl)); out["worthy_clones"] = int(cl.worthy.sum())
    out["true_accepts"] = int(((cl.decision == "accept") & cl.worthy).sum()); out["false_accepts"] = int(((cl.decision == "accept") & ~cl.worthy).sum())
    out["resequence_rate"] = float((cl.decision == "re-sequence").mean())
    out["decision_counts"] = cl.decision.value_counts().to_dict()
    cl[(cl.decision == "accept") & ~cl.worthy][["i", "project", "scenario", "truth_goal"]].to_csv(a.prefix + "_false_accepts.csv", index=False)
    pd.crosstab(cl.decision, cl.worthy).to_csv(a.prefix + "_decision_matrix.csv")
    cl[(cl.decision != "accept") & cl.worthy][["i", "project", "scenario", "truth_goal", "decision", "tier"]].to_csv(a.prefix + "_missed_accepts.csv", index=False)
    json.dump(out, open(a.prefix + "_extras.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
    print(tiers.to_string())


if __name__ == "__main__":
    main()
