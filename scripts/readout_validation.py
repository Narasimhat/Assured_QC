#!/usr/bin/env python
"""Independent readout of Sanger traces for validating a quantification tool without running it.

For each control/edited pair it measures, directly on the chromatograms:
  * SNP / ssODN jobs: the share of the donor base against the wild-type base at every position where the donor differs from the
    control read (dye-scaled peak heights at the called peak).  Edits that are present on all alleles read ~1, absent ~0, one allele ~0.5.
  * all jobs: the single-sequence purity (top dye share) downstream of the first cut relative to the same read upstream, and the
    in-phase base-call identity to the control downstream of the cut.  A single allele that differs from the wild type reads
    purity ratio ~1 with low identity; a wild-type or a clean mixture of alleles reads differently.
Usage: python readout_validation.py jobs.json out.csv
"""
import sys, json, numpy as np, pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
BASES = "ACGT"
rc = lambda s: str(Seq(s).reverse_complement())
al = Align.PairwiseAligner(); al.mode = "local"; al.match_score = 2; al.mismatch_score = -2; al.open_gap_score = -6; al.extend_gap_score = -2

def load(path):
    rec = SeqIO.read(path, "abi"); raw = rec.annotations["abif_raw"]
    order = raw["FWO_1"].decode(); ch = {b: np.array(raw[f"DATA{9 + i}"], dtype=float) for i, b in enumerate(order)}
    ploc = np.array(raw["PLOC2"]); seq = str(rec.seq); L = min(len(seq), len(ploc))
    h = np.array([[ch[b][max(0, p - 1):p + 2].max() for b in BASES] for p in ploc[:L]])
    sc = []
    for b in range(4):
        sel = [i for i in range(L) if seq[i] == BASES[b] and h[i].sum() > 200]
        sc.append(np.median(h[sel, b]) if sel else 1.0)
    sc = np.array(sc) / np.mean(sc); hn = h / sc; tot = hn.sum(1)
    return dict(seq=seq[:L], comp=hn / np.maximum(tot[:, None], 1), tot=tot)

def locate(seq, query):
    """Return (score, strand, alignment) of the best local hit of query or its reverse complement."""
    query = query.upper(); best = None
    for strand, q in (("+", query), ("-", rc(query))):
        a = al.align(seq, q)[0]
        if best is None or a.score > best[0]: best = (a.score, strand, a)
    return best

def offset_between(ctrl, edit, lo, hi):
    """Offset (edit index - control index) at the control stretch lo..hi, taken from a local alignment of the two base-call strings,
    and the identity of that stretch (so a drifting early-read register does not matter)."""
    a = al.align(ctrl["seq"], edit["seq"])[0]; m = {}
    for (t0, t1), (q0, q1) in zip(*a.aligned):
        for k in range(t1 - t0): m[t0 + k] = q0 + k
    ps = [p for p in range(max(0, lo), hi) if p in m]
    if len(ps) < 8: return 0, 0.0
    offs = [m[p] - p for p in ps]; off = int(np.median(offs))
    return off, float(np.mean([ctrl["seq"][p] == edit["seq"][p + off] for p in ps if 0 <= p + off < len(edit["seq"])]))

def purity(r, lo, hi):
    j = [i for i in range(max(0, lo), min(len(r["seq"]), hi)) if r["tot"][i] > 100]
    return (float(np.mean(r["comp"][j].max(1))) if len(j) >= 20 else np.nan), len(j)

def readout(job):
    out = {}
    c, e = load(job["control"]), load(job["edited"])
    # cut positions on the control read for each guide (3 bp from the PAM-proximal end of the protospacer)
    cuts = []
    for g in job["guides"]:
        sc, strand, a = locate(c["seq"], g)
        (ts, qs) = a.aligned
        if sc < 30: continue
        s0, s1 = int(ts[0][0]), int(ts[-1][1])
        cuts.append((s0 + 17, strand, sc) if strand == "+" else (s1 - 17, strand, sc))
    out["n_cuts_found"] = len(cuts)
    if not cuts: return out
    # reads run in one direction: "downstream" of the first cut = increasing index
    cut1 = min(x[0] for x in cuts); out["cut1_read_pos"] = cut1
    off, idn = offset_between(c, e, cut1 - 70, cut1 - 20); out["offset"] = off; out["upstream_identity"] = round(idn, 3)
    pu, nu = purity(e, cut1 - 110 + off, cut1 - 10 + off); pd_, nd = purity(e, cut1 + 10 + off, cut1 + 110 + off)
    out["purity_up"] = pu; out["purity_down"] = pd_; out["purity_ratio"] = pd_ / pu if pu and pu == pu and pd_ == pd_ else np.nan
    idx = [i for i in range(cut1 + 10, cut1 + 110) if 0 <= i < len(c["seq"]) and 0 <= i + off < len(e["seq"]) and c["tot"][i] > 100 and e["tot"][i + off] > 100]
    out["inphase_identity_down"] = float(np.mean([c["seq"][i] == e["seq"][i + off] for i in idx])) if len(idx) >= 30 else np.nan
    if job["donor"]:
        sc, strand, a = locate(c["seq"], job["donor"]); donor = job["donor"].upper() if strand == "+" else rc(job["donor"].upper())
        (ts, qs) = a.aligned; fr = []
        for (t0, t1), (q0, q1) in zip(ts, qs):
            for k in range(t1 - t0):
                p = t0 + k; db = donor[q0 + k]; cb = c["seq"][p]
                if db != cb and db in BASES and cb in BASES and c["comp"][p][BASES.index(cb)] > 0.6 and 0 <= p + off < len(e["seq"]):
                    cp = e["comp"][p + off]; den = cp[BASES.index(db)] + cp[BASES.index(cb)]
                    if den > 0.6 and e["tot"][p + off] > 100: fr.append((p, db, cb, float(cp[BASES.index(db)] / den)))
        out["donor_sites"] = json.dumps([(int(p), db, cb, round(f, 2)) for p, db, cb, f in fr]); out["n_donor_sites"] = len(fr)
        if fr: out["donor_frac_min"] = min(x[3] for x in fr); out["donor_frac_max"] = max(x[3] for x in fr); out["donor_frac_mean"] = float(np.mean([x[3] for x in fr]))
    return out

if __name__ == "__main__":
    jobs = json.load(open(sys.argv[1])); rows = []
    for i, j in enumerate(jobs):
        try: r = readout(j)
        except Exception as ex: r = {"error": str(ex)}
        r["idx"] = i; rows.append(r)
    pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
