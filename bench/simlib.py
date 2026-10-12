"""In-silico truth generator for Sanger edit-QC benchmarks.

Edited traces are built from REAL control traces (real peak shapes, heights, context effects, noise and run quality):
each allele is rendered by copying, base by base, the sample windows of the control trace around the control peak that
carries the same sequence context (deleted bases are skipped, so downstream peaks are the control's own peaks moved into
place); bases that do not exist in the control (substitutions, insertions) borrow a nearby control peak of that base.
The mixture is the fraction-weighted sum of the allele renderings, followed by optional baseline noise, per-channel dye
bias and an overall amplitude change. The truth (allele sizes, kinds, fractions) is exact by construction.

Limits (stated, not hidden): PCR length bias is not modelled (truth = signal fraction); the reverse read is the forward
trace reversed with fresh noise, not an independent run; sequence context effects of novel junctions are borrowed from
the control. These traces test the engine; they do not replace the wet-lab mixtures.
"""
import os, json, struct, numpy as np
from Bio import SeqIO, Align

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s[::-1].translate(COMP)

def _entry(name, num, etype, esize, n, data): return dict(name=name, num=num, etype=etype, esize=esize, n=n, data=data)
def build_abif(entries):
    out = bytearray(b"ABIF" + struct.pack(">H", 101) + b"\0" * 28); out += b"\0" * (128 - len(out)); placed = []
    for e in entries:
        data = e["data"]
        if len(data) <= 4: off = int.from_bytes(data.ljust(4, b"\0"), "big")
        else:
            off = len(out); out += data
            if len(out) % 2: out += b"\0"
        placed.append((e, off))
    dir_off = len(out)
    for e, off in placed: out += struct.pack(">4siHHiiii", e["name"].encode(), e["num"], e["etype"], e["esize"], e["n"], len(e["data"]), off, 0)
    out[6:34] = struct.pack(">4siHHiiii", b"tdir", 1, 1023, 28, len(entries), 28 * len(entries), dir_off, 0)
    return bytes(out)

def write_ab1(path, chans, calls, ploc, qual, name="sim"):
    """chans: (4, n) heights in ACGT order."""
    FWO = "GATC"; order = ["ACGT".index(b) for b in FWO]
    data = [np.clip(chans[o], -32000, 32000).round().astype(">i2").tobytes() for o in order]
    L = len(calls); cb = calls.encode(); qb = np.asarray(qual, dtype=np.uint8).tobytes(); pl = np.asarray(ploc).astype(">i2").tobytes()
    ents = [_entry("DATA", 9 + k, 4, 2, chans.shape[1], data[k]) for k in range(4)]
    ents += [_entry("FWO_", 1, 2, 1, 4, FWO.encode()), _entry("PBAS", 1, 2, 1, L, cb), _entry("PBAS", 2, 2, 1, L, cb),
             _entry("PCON", 1, 2, 1, L, qb), _entry("PCON", 2, 2, 1, L, qb), _entry("PLOC", 1, 4, 2, L, pl), _entry("PLOC", 2, 4, 2, L, pl),
             _entry("SMPL", 1, 18, 1, len(name) + 1, bytes([len(name)]) + name.encode())]
    open(path, "wb").write(build_abif(sorted(ents, key=lambda e: (e["name"], e["num"]))))

class Control:
    def __init__(self, path=None, arrays=None):
        if arrays is not None:
            self.path = None; self.ch, self.calls, self.ploc, self.pcon = arrays
        else:
            rec = SeqIO.read(path, "abi"); ar = rec.annotations["abif_raw"]; self.path = path
            fwo = ar["FWO_1"]; fwo = fwo.decode() if isinstance(fwo, bytes) else fwo
            ch = {fwo[k]: np.array(ar[f"DATA{9 + k}"], dtype=float) for k in range(4)}
            self.ch = np.stack([ch[b] for b in "ACGT"]); calls = ar["PBAS2"]; self.calls = (calls.decode() if isinstance(calls, bytes) else calls).upper()
            self.ploc = np.array(ar["PLOC2"], dtype=int); pc = ar["PCON2"]; self.pcon = np.frombuffer(pc, dtype=np.uint8).astype(int) if isinstance(pc, bytes) else np.array(pc, dtype=int)
        self.L, self.n = len(self.calls), self.ch.shape[1]
        top = np.array([self.ch[:, max(0, p - 1):p + 2].max(axis=1).sum() for p in self.ploc])
        from scipy.ndimage import median_filter
        self.env = np.maximum(median_filter(top, size=41, mode="nearest"), 1.0)
        head = top[:min(self.L, 700)]; self.typical = float(np.median(head[head > 100])) if (head > 100).any() else float(np.median(top))
        strong = np.where(self.env >= 0.4 * np.percentile(self.env[:min(self.L, 700)], 90))[0]
        # usable_end: end of the first run of strong signal (the read's good region); sources beyond it are noise, not peaks
        gaps = np.where(np.diff(strong) > 20)[0]; self.usable_end = int(strong[gaps[0]] if len(gaps) else strong[-1]) if len(strong) else self.L - 1
        self.lo = np.empty(self.L, dtype=int); self.hi = np.empty(self.L, dtype=int)
        mid = (self.ploc[:-1] + self.ploc[1:]) // 2
        self.lo[1:], self.hi[:-1] = mid, mid
        self.lo[0] = max(0, self.ploc[0] - (self.ploc[1] - self.ploc[0]) // 2); self.hi[-1] = min(self.n, self.ploc[-1] + (self.ploc[-1] - self.ploc[-2]) // 2)

    def truncated(self, length):
        """The first `length` bases: the amplicon as a shorter product (a read cannot be longer than the product it reads)."""
        hi = int(self.hi[length - 1]); return Control(arrays=(self.ch[:, :hi].copy(), self.calls[:length], self.ploc[:length].copy(), self.pcon[:length].copy()))

    def reversed(self):
        """The same template read from the other primer: channels swapped (A<->T, C<->G) and the trace run backwards."""
        n = self.n; return Control(arrays=(self.ch[[3, 2, 1, 0]][:, ::-1].copy(), rc(self.calls), (n - 1 - self.ploc)[::-1].copy(), self.pcon[::-1].copy()))

    def find_guide(self, g, mm=1):
        best = None
        for strand, s in (("+", g), ("-", rc(g))):
            for i in range(self.L - len(s) + 1):
                d = sum(a != b for a, b in zip(self.calls[i:i + len(s)], s))
                if d <= mm and (best is None or d < best[0]): best = (d, i, strand)
        if best is None: raise ValueError("guide not in control read")
        return best[1] + 17 if best[2] == "+" else best[1] + 3   # index of the first base after the cut

    def donor_variants(self, donor):
        al = Align.PairwiseAligner(); al.mode = "local"; al.match_score = 2; al.mismatch_score = -3; al.open_gap_score = -6; al.extend_gap_score = -2
        best = None
        for d in (donor, rc(donor)):
            a = al.align(self.calls, d)[0]
            if best is None or a.score > best[0].score: best = (a, d)
        a, d = best; vs = []
        for (r0, r1), (d0, d1) in zip(*a.aligned):
            for k in range(r1 - r0):
                if self.calls[r0 + k] != d[d0 + k]: vs.append((int(r0 + k), d[d0 + k]))
        return vs

def build_allele(L, ops):
    """ops: ('sub', pos, base) | ('del', pos, n) | ('ins', pos, seq); positions are control indices; applied right to left."""
    items = list(range(L))
    for op in sorted(ops, key=lambda o: -o[1]):
        k = items.index(op[1])
        if op[0] == "sub": items[k] = (op[2], op[1])
        elif op[0] == "del": del items[k:k + op[2]]
        elif op[0] == "ins": items[k:k] = [(b, op[1]) for b in op[2]]
    return items

def allele_seq(ctl, items): return "".join(ctl.calls[i] if isinstance(i, (int, np.integer)) else i[0] for i in items)

def _borrow(ctl, base, anchor, cache):
    key = (base, anchor)
    if key not in cache:
        cand = [m for m in range(max(0, anchor - 25), min(ctl.L, anchor + 26)) if ctl.calls[m] == base and ctl.pcon[m] >= 25 and abs(m - anchor) > 1]
        if not cand: cand = [m for m in range(ctl.L) if ctl.calls[m] == base]
        cache[key] = min(cand, key=lambda m: abs(m - anchor))
    return cache[key]

def render(ctl, alleles, rng, noise=0.0, dye=(1, 1, 1, 1), amp=1.0, jitter=0.0):
    """alleles: list of (items, fraction). Returns channels (4, n_out), calls, ploc, qual."""
    def usable(items):
        for j, it in enumerate(items):
            if isinstance(it, (int, np.integer)) and it > ctl.usable_end: return j
        return len(items)
    Lout = min(ctl.L, min(min(len(it), usable(it)) for it, _ in alleles)); cache = {}
    n_out = int(ctl.hi[Lout - 1]); out = np.zeros((4, n_out)); out[:, :ctl.lo[0]] = ctl.ch[:, :ctl.lo[0]]
    comp = np.zeros((Lout, 4))
    for items, f in alleles:
        for j in range(Lout):
            item = items[j]; lo, hi = ctl.lo[j], ctl.hi[j]
            if isinstance(item, (int, np.integer)): m, base = int(item), ctl.calls[int(item)]
            else: base = item[0]; m = _borrow(ctl, base, item[1], cache)
            if m == j: out[:, lo:hi] += f * ctl.ch[:, lo:hi]
            else:
                scale = float(np.clip(ctl.env[j] / ctl.env[m], 0.4, 2.5)); idx = np.clip(ctl.ploc[m] + np.arange(lo - ctl.ploc[j], hi - ctl.ploc[j]), 0, ctl.n - 1)
                out[:, lo:hi] += f * scale * ctl.ch[:, idx]
            if base in "ACGT": comp[j, "ACGT".index(base)] += f
    if jitter > 0:   # run-to-run variation: per-peak height jitter (per channel too) and a smooth gain drift along the read
        drift = np.exp(np.interp(np.arange(Lout), np.linspace(0, Lout - 1, 8), rng.normal(0, 1.5 * jitter, 8)))
        for j in range(Lout):
            out[:, ctl.lo[j]:ctl.hi[j]] *= (drift[j] * np.exp(rng.normal(0, jitter)) * np.exp(rng.normal(0, 0.5 * jitter, 4)))[:, None]
    if noise > 0: out += rng.normal(0, noise * ctl.typical, out.shape)
    out *= (np.asarray(dye, dtype=float) * amp)[:, None]
    calls = "".join("ACGT"[int(np.argmax(r))] for r in comp)
    main = comp.max(axis=1); qual = np.clip(3 * main * 100 - 240, 0, 60).astype(int)
    same = np.array([all(isinstance(it[j], (int, np.integer)) and int(it[j]) == j for it, _ in alleles) for j in range(Lout)])
    qual = np.where(same, ctl.pcon[:Lout], np.minimum(qual, ctl.pcon[:Lout] if False else 60))
    return out, calls, ctl.ploc[:Lout].copy(), qual

def reverse_read(chans, calls, ploc, qual, rng, noise=0.01, typical=1500.0):
    swap = [3, 2, 1, 0]                       # A<->T, C<->G
    rev = chans[swap][:, ::-1].copy(); n = rev.shape[1]
    rev += rng.normal(0, noise * typical, rev.shape)
    return rev, rc(calls), (n - 1 - ploc)[::-1].copy(), np.asarray(qual)[::-1].copy()
