"""Extra synthetic fixtures for events outside the exhaustive small-indel library: a large deletion, a 6-base insertion of unknown
sequence, a wild-type replicate, and an ambiguous phase (wild type + edit-with-indel). Built from fixtures/sim/truth.json (references,
guides, cuts and the existing controls); all traces are generated (sim_ab1.py), no third-party data.
Usage: python make_extended_fixtures.py fixtures/sim
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from sim_ab1 import *

OUT = sys.argv[1]; truth = json.load(open(os.path.join(OUT, "truth.json")))
def fit(seq, L, seed):
    rng = np.random.default_rng(seed); return seq[:L] if len(seq) >= L else seq + "".join(rng.choice(list("ACGT"), L - len(seq)))
def a_del(ref, c, n): k = n // 2; return ref[:c - k] + ref[c - k + n:]
def a_ins(ref, c, s): return ref[:c] + s + ref[c:]
def write_edited(path, ref, alleles, seed, noise_frac=0.01, noise=0.012):
    L = len(ref); seqs = [(lab, fit(s, L, seed + k), f) for k, (lab, s, f) in enumerate(alleles)]
    assert abs(sum(f for _, _, f in seqs) - 1) < 1e-9
    sig, comp = signals_from_alleles([(s, f) for _, s, f in seqs], amplitude_profile(L, seed=1), noise_frac=noise_frac, seed=seed)
    write_ab1(path, "".join(BASES[int(np.argmax(r))] for r in comp), sig, seed=seed, noise=noise)

new = {}
nr = truth["nr2f2_ko_clone_het_wt_m1"]; ref, c1 = nr["reference"], nr["cuts"]["g1"]
ap = truth["apoe_snp_clone_hom"]; aref, ac1 = ap["reference"], ap["cuts"]["g1"]; hdr = ap["alleles"][0]["sequence"]
cases = {
    "ext_nr2f2_large_del60_het": (ref, nr, [("wt", ref, .5), ("del-60", a_del(ref, c1, 60), .5)], "ko_clone_large_del"),
    "ext_nr2f2_ins6_het": (ref, nr, [("wt", ref, .5), ("ins+6", a_ins(ref, c1, "GATTAC"), .5)], "ko_clone_ins6"),
    "ext_nr2f2_abe_p8_40": (ref, nr, [("wt", ref, .6), ("A>G p8", (lambda g: ref[:ref.find(g) + 7] + "G" + ref[ref.find(g) + 8:])(nr["guides"][0]), .4)], "base_edit_abe_p8_40pct"),
    "ext_nr2f2_wt_replicate": (ref, nr, [("wt", ref, 1.0)], "wild_type"),
    "ext_apoe_flanked": (aref, ap, [("wt", aref, .5), ("hdr+ins1", a_ins(hdr, ac1, "A"), .5)], "hdr_with_indel_wt_markers_flank_the_indel"),
    # one SNP 15 bases upstream of the cut, one guide (as in a single-SNP ssODN design): the insertion is downstream of every marker, so the phase cannot be read
    "ext_apoe_phase_ambiguous": (aref, ap, [("wt", aref, .5), ("snp+ins1", a_ins(aref[:ac1 - 15] + ("C" if aref[ac1 - 15] != "C" else "G") + aref[ac1 - 14:], ac1, "A"), .5)], "snp_with_indel_wt_phase_ambiguous"),
}
for k, (name, (r, base, alleles, kind)) in enumerate(cases.items()):
    write_edited(os.path.join(OUT, name + "_edited.ab1"), r, alleles, 500 + k, *((0.03, 0.03) if name == "ext_apoe_phase_ambiguous" else ()))
    donor = None; guides = base["guides"]
    if name == "ext_apoe_phase_ambiguous":
        c = ac1; snp = ("C" if aref[c - 15] != "C" else "G"); donor = aref[c - 60:c - 15] + snp + aref[c - 14:c + 60]; guides = base["guides"][:1]
    new[name] = dict(control=base["control"], edited=name + "_edited.ab1", guides=guides, donor=donor or base["donor"], kind=kind, reference=r, cuts=base["cuts"],
                     alleles=[dict(label=l, fraction=f, indel=len(q) - len(r)) for l, q, f in alleles])
json.dump(new, open(os.path.join(OUT, "truth_extended.json"), "w"), indent=1); print(sorted(new))
