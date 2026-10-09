"""Build the synthetic QC fixtures: control and edited .ab1 files with known allele mixtures.
Usage: python make_qc_fixtures.py <reference fixtures dir> <APOE design summary.json> <NR2F2 small-tag design json> <outdir> [<NR2F2 block design spec json>]
All traces are generated here (sim_ab1.py); no third-party trace data is used."""
import sys, os, json, hashlib, numpy as np
from Bio import SeqIO
from sim_ab1 import *

FXD, SNPJ, TAGJ, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)

def rcseq(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))
def cut_index(read, guide):
    i = read.find(guide)
    if i >= 0: return i + 17
    j = read.find(rcseq(guide))
    if j >= 0: return j + 3
    raise ValueError("guide not in read")
def fit(seq, L, seed):
    rng = np.random.default_rng(seed)
    return seq[:L] if len(seq) >= L else seq + "".join(rng.choice(list("ACGT"), L - len(seq)))
def a_del(ref, c, n): k = n // 2; return ref[:c-k] + ref[c-k+n:]
def a_ins(ref, c, s): return ref[:c] + s + ref[c:]
def a_hdr(ref, donor, la=30, ra=20):
    iL, iR = ref.find(donor[:la]), ref.find(donor[-ra:]); assert iL >= 0 and iR >= 0
    return ref[:iL] + donor + ref[iR+ra:]

def write_control(path, ref, seed):
    sig, _ = signals_from_alleles([(ref, 1.0)], amplitude_profile(len(ref), seed=1), seed=seed)
    write_ab1(path, ref, sig, seed=seed)
def write_edited(path, ref, alleles, seed):
    L = len(ref); seqs = [(lab, fit(s, L, seed + k), f) for k, (lab, s, f) in enumerate(alleles)]
    assert abs(sum(f for _, _, f in seqs) - 1) < 1e-9
    sig, comp = signals_from_alleles([(s, f) for _, s, f in seqs], amplitude_profile(L, seed=1), seed=seed)
    write_ab1(path, "".join(BASES[int(np.argmax(r))] for r in comp), sig, seed=seed)
    return seqs

nr = str(SeqIO.read(os.path.join(FXD, "nr2f2-ng016753.gb"), "genbank").seq).upper()
ap = str(SeqIO.read(os.path.join(FXD, "apoe-r154s.gb"), "genbank").seq).upper()
tag = json.load(open(TAGJ)); snp = json.load(open(SNPJ))["two"]
g1, g2 = "CAGTTTTAACTGGCCGTATA", "AATAAATAAATAAAATAAGA"
ga1, ga2 = "ACACTGCCAGGCGCTTCTGC", "CCCCGGCCTGGTACACTGCC"
ref = nr[nr.find("ACACACCTCATGTGACCC"):][:760]
aref = ap[ap.find(snp["primers"][0]["s"]):][:760]
c1, c2 = cut_index(ref, g1), cut_index(ref, g2); ac1, ac2 = cut_index(aref, ga1), cut_index(aref, ga2)
od1, wo1 = snp["donors"][0]["od"], snp["donors"][0]["wo"]
hdr1 = a_hdr(aref, od1, 40, 18); silent_only = a_hdr(aref, od1[:95] + wo1[95] + od1[96:], 40, 18)
ha = a_hdr(ref, tag["donor"])
S = {}
def add(name, r, ctl, alleles, guides, donor=None, kind=""): S[name] = dict(ref=r, ctl=ctl, alleles=alleles, guides=guides, donor=donor, kind=kind)
add("nr2f2_ko_pool_g1", ref, "nr2f2_control.ab1", [("wt",ref,.25),("del-1",a_del(ref,c1,1),.20),("ins+1",a_ins(ref,c1,ref[c1-1]),.25),("del-7",a_del(ref,c1,7),.08),("del-12",a_del(ref,c1,12),.12),("del-3",a_del(ref,c1,3),.10)], [g1], kind="ko_pool")
add("nr2f2_ko_clone_het_wt_m1", ref, "nr2f2_control.ab1", [("wt",ref,.5),("del-1",a_del(ref,c1,1),.5)], [g1], kind="ko_clone")
add("nr2f2_ko_clone_compound", ref, "nr2f2_control.ab1", [("ins+1",a_ins(ref,c1,ref[c1-1]),.5),("del-7",a_del(ref,c1,7),.5)], [g1], kind="ko_clone")
add("nr2f2_ko_clone_hom_m1", ref, "nr2f2_control.ab1", [("del-1",a_del(ref,c1,1),1.0)], [g1], kind="ko_clone")
add("nr2f2_del_two_guides", ref, "nr2f2_control.ab1", [("wt",ref,.30),("del_between_cuts",ref[:c1]+ref[c2:],.25),("g1:del-1",a_del(ref,c1,1),.15),("g2:ins+1",a_ins(ref,c2,ref[c2-1]),.15),("g2:del-2",a_del(ref,c2,2),.15)], [g1,g2], kind="deletion")
add("nr2f2_ha_hdr_pool", ref, "nr2f2_control.ab1", [("hdr",ha,.30),("wt",ref,.30),("g1:del-1",a_del(ref,c1,1),.15),("g1:ins+1",a_ins(ref,c1,ref[c1-1]),.10),("g2:del-2",a_del(ref,c2,2),.15)], [g1,g2], tag["donor"], "hdr_pool")
add("nr2f2_ha_hdr_clone_het", ref, "nr2f2_control.ab1", [("hdr",ha,.5),("wt",ref,.5)], [g1,g2], tag["donor"], "hdr_clone")
add("nr2f2_ha_hdr_clone_hom", ref, "nr2f2_control.ab1", [("hdr",ha,1.0)], [g1,g2], tag["donor"], "hdr_clone")
add("apoe_snp_hdr_pool", aref, "apoe_control.ab1", [("hdr",hdr1,.35),("silent_only",silent_only,.10),("wt",aref,.20),("g1:del-1",a_del(aref,ac1,1),.15),("g1:ins+1",a_ins(aref,ac1,aref[ac1-1]),.10),("g2:del-4",a_del(aref,ac2,4),.10)], [ga1,ga2], od1, "hdr_pool")
add("apoe_snp_clone_het", aref, "apoe_control.ab1", [("hdr",hdr1,.5),("wt",aref,.5)], [ga1,ga2], od1, "hdr_clone")
add("apoe_snp_clone_hom", aref, "apoe_control.ab1", [("hdr",hdr1,1.0)], [ga1,ga2], od1, "hdr_clone")
add("apoe_snp_clone_partial_het", aref, "apoe_control.ab1", [("hdr",hdr1,.5),("silent_only",silent_only,.5)], [ga1,ga2], od1, "hdr_clone")

# reverse-primer read of the same locus: the read starts at out_R and runs towards the stop codon
iR = nr.find(rcseq("CCATTCTGCTAATTGTCTCCC")) + len("CCATTCTGCTAATTGTCTCCC")
ref_rev = rcseq(nr[iR - 760: iR]); rc1 = cut_index(ref_rev, g1)
add("nr2f2_rev_ko_pool_g1", ref_rev, "nr2f2_rev_control.ab1", [("wt",ref_rev,.25),("del-1",a_del(ref_rev,rc1,1),.20),("ins+1",a_ins(ref_rev,rc1,ref_rev[rc1-1]),.25),("del-7",a_del(ref_rev,rc1,7),.08),("del-12",a_del(ref_rev,rc1,12),.12),("del-3",a_del(ref_rev,rc1,3),.10)], [g1], kind="ko_pool_reverse_read")

write_control(os.path.join(OUT, "nr2f2_control.ab1"), ref, 1)
write_control(os.path.join(OUT, "nr2f2_rev_control.ab1"), ref_rev, 3)
write_control(os.path.join(OUT, "apoe_control.ab1"), aref, 2)
truth = {}
for k, (name, s) in enumerate(S.items()):
    write_edited(os.path.join(OUT, name + "_edited.ab1"), s["ref"], s["alleles"], 100 + k)
    truth[name] = dict(control=s["ctl"], edited=name + "_edited.ab1", guides=s["guides"], donor=s["donor"], kind=s["kind"],
                       reference=s["ref"], cuts=dict(g1=int(c1), g2=int(c2)) if s["ref"] is ref else (dict(g1=int(rc1)) if s["ref"] is ref_rev else dict(g1=int(ac1), g2=int(ac2))),
                       alleles=[dict(label=l, fraction=f, indel=len(q) - len(s["ref"]), sequence=fit(q, len(s["ref"]), 0)) for l, q, f in s["alleles"]])
json.dump(truth, open(os.path.join(OUT, "truth.json"), "w"), indent=1)
print({n: hashlib.sha256(open(os.path.join(OUT, n), "rb").read()).hexdigest()[:10] for n in sorted(os.listdir(OUT)) if n.endswith(".ab1")})


# ---- artefact and noise stress samples (APOE and NR2F2 loci)
def artifact_profile(L, seed, centre=None, n=45, lo=0.25, hi=0.45):
    # shoulders and uneven peaks, concentrated in the stretch the analysis reads (from 25 bases before the cut to 120 after)
    rng = np.random.default_rng(seed)
    pos = rng.choice(np.arange(centre - 25, centre + 120), n, replace=False)
    return {int(p): (int(rng.integers(0, 4)), float(rng.uniform(lo, hi))) for p in pos}
def write_with(path, ref_, alleles, seed, artifact=None, noise_frac=0.01, noise=0.012):
    L = len(ref_); seqs = [(lab, fit(s, L, seed + k), f) for k, (lab, s, f) in enumerate(alleles)]
    sig, comp = signals_from_alleles([(s, f) for _, s, f in seqs], amplitude_profile(L, seed=1), noise_frac=noise_frac, seed=seed, artifact=artifact, context_ref=ref_)
    write_ab1(path, "".join(BASES[int(np.argmax(r))] for r in comp), sig, seed=seed, noise=noise)
    return seqs
art = artifact_profile(len(ref), 7, centre=c1)
write_with(os.path.join(OUT, "nr2f2_artifact_control.ab1"), ref, [("wt", ref, 1.0)], 41, artifact=art)
extra = {}
for k, (name, alleles) in enumerate([
    ("nr2f2_artifact_null", [("wt", ref, 1.0)]),
    ("nr2f2_artifact_ko_het", [("wt", ref, .5), ("del-1", a_del(ref, c1, 1), .5)]),
]):
    write_with(os.path.join(OUT, name + "_edited.ab1"), ref, alleles, 410 + k, artifact=art)
    extra[name] = dict(control="nr2f2_artifact_control.ab1", edited=name + "_edited.ab1", guides=[g1], donor=None, kind="artifact", reference=ref, cuts=dict(g1=int(c1), g2=int(c2)),
                       alleles=[dict(label=l, fraction=f, indel=len(q) - len(ref), sequence=fit(q, len(ref), 0)) for l, q, f in alleles])
write_with(os.path.join(OUT, "apoe_noisy_control.ab1"), aref, [("wt", aref, 1.0)], 51, noise_frac=0.02, noise=0.025)
for k, (name, alleles) in enumerate([
    ("apoe_noisy_snp_het", [("hdr", hdr1, .5), ("wt", aref, .5)]),
    ("apoe_noisy_snp_hom", [("hdr", hdr1, 1.0)]),
]):
    write_with(os.path.join(OUT, name + "_edited.ab1"), aref, alleles, 420 + k, noise_frac=0.02, noise=0.025)
    extra[name] = dict(control="apoe_noisy_control.ab1", edited=name + "_edited.ab1", guides=[ga1, ga2], donor=od1, kind="noisy_hdr_clone", reference=aref, cuts=dict(g1=int(ac1), g2=int(ac2)),
                       alleles=[dict(label=l, fraction=f, indel=len(q) - len(aref), sequence=fit(q, len(aref), 0)) for l, q, f in alleles])
truth.update(extra)
json.dump(truth, open(os.path.join(OUT, "truth.json"), "w"), indent=1)

# ---- junction reads (no control): the SD40-V5 block knock-in allele read across its 5' and 3' junctions
BLOCK_SPEC = sys.argv[5] if len(sys.argv) > 5 else ""
spec_blk = json.load(open(BLOCK_SPEC)) if BLOCK_SPEC and os.path.exists(BLOCK_SPEC) else None
if spec_blk is not None:
    W = spec_blk["reference"]["sequence"]; d0 = spec_blk["donors"][0]
    allele = W[:d0["refStart"]] + d0["sequence"] + W[d0["refEnd"]:]
    ins0 = d0["refStart"] + d0["arm5"]                       # first base of the cassette in allele coordinates
    fw = spec_blk["primers"][0]; rv = spec_blk["primers"][1]
    read5 = allele[fw["start"]: fw["start"] + 760]
    rev_end = rv["end"] + (len(allele) - len(W))               # reverse primer end in allele coordinates
    read3 = rcseq(allele[rev_end - 760: rev_end])
    def mutate(seq, i, new): return seq[:i] + new + seq[i+1:]
    jr = {}
    jr["junction5_ok"] = read5
    off = ins0 - fw["start"]                                    # cassette start inside read5
    jr["junction5_del_in_insert"] = read5[:off + 100] + read5[off + 101:] + "A"
    arm_i = off - 120
    jr["junction5_arm_snp"] = mutate(read5, arm_i, "A" if read5[arm_i] != "A" else "C")
    mk = [m for m in spec_blk["markers"]][0]["pos"] - fw["start"]   # first blocking change inside read5
    jr["junction5_blocking_missing"] = mutate(read5, mk, W[spec_blk["markers"][0]["pos"]])
    jr["junction3_ok"] = read3
    jr_truth = {}
    for k, (name, seq) in enumerate(jr.items()):
        L = len(seq)
        sig, comp = signals_from_alleles([(seq, 1.0)], amplitude_profile(L, seed=1), seed=300 + k)
        write_ab1(os.path.join(OUT, name + ".ab1"), seq, sig, seed=300 + k)
        jr_truth[name] = dict(file=name + ".ab1", read=seq)
    jr_truth["_meta"] = dict(allele_length=len(allele), insert_start_in_allele=ins0, read5_cassette_offset=off, arm_snp_index=arm_i, blocking_marker_index=mk)
    json.dump(jr_truth, open(os.path.join(OUT, "junction_truth.json"), "w"), indent=1)
