
"""Synthetic Sanger trace generator (own code, no third-party trace data).
Produces ABIF files carrying the tags the common readers need: DATA9-12, FWO_1, PBAS1/2, PLOC1/2, PCON1/2, SMPL1.
The edited-sample trace follows the linear mixture model: at each base position the four channel heights are
the total signal times the sum of the allele fractions carrying that base.
"""
import struct, numpy as np

BASES = "ACGT"
FWO = "GATC"
CH_SCALE = {"A": 1.0, "C": 0.85, "G": 0.7, "T": 0.95}   # channel-specific height factors
SPACING = 11
SIGMA = 2.3

def _entry(name, num, etype, esize, n, data):
    return dict(name=name, num=num, etype=etype, esize=esize, n=n, data=data)

def build_abif(entries):
    out = bytearray(b"ABIF" + struct.pack(">H", 101) + b"\0" * 28)
    out += b"\0" * (128 - len(out))
    placed = []
    for e in entries:
        data = e["data"]
        if len(data) <= 4:
            off = int.from_bytes(data.ljust(4, b"\0"), "big")
        else:
            off = len(out); out += data
            if len(out) % 2: out += b"\0"
        placed.append((e, off))
    dir_off = len(out)
    for e, off in placed:
        out += struct.pack(">4siHHiiii", e["name"].encode(), e["num"], e["etype"], e["esize"], e["n"], len(e["data"]), off, 0)
    out[6:34] = struct.pack(">4siHHiiii", b"tdir", 1, 1023, 28, len(entries), 28 * len(entries), dir_off, 0)
    return bytes(out)

def make_trace(base_signals, seed=0, noise=0.012, background=12.0, quality_floor=None):
    """base_signals: (L, 4) array of channel heights in ACGT order per base position. Returns channel arrays + ploc."""
    rng = np.random.default_rng(seed)
    L = base_signals.shape[0]
    n = SPACING * (L + 4)
    x = np.arange(n)
    chans = np.zeros((4, n))
    ploc = np.array([SPACING * (j + 2) + int(round(rng.normal(0, 0.4))) for j in range(L)])
    for j in range(L):
        lo, hi = max(0, ploc[j] - 5 * SPACING // 2), min(n, ploc[j] + 5 * SPACING // 2)
        g = np.exp(-0.5 * ((x[lo:hi] - ploc[j]) / SIGMA) ** 2)
        for c in range(4):
            chans[c, lo:hi] += base_signals[j, c] * g
    chans += rng.normal(0, background, chans.shape).clip(min=0)
    chans *= (1 + rng.normal(0, noise, chans.shape))
    return chans, ploc

def amplitude_profile(L, seed=0, start_poor=45, end_decay=80):
    rng = np.random.default_rng(seed + 101)
    amp = 1500 * np.exp(rng.normal(0, 0.12, L))
    j = np.arange(L)
    amp *= np.where(j < start_poor, 0.25 + 0.75 * j / start_poor, 1.0)
    amp *= np.where(j > L - end_decay, np.maximum(0.15, 1 - (j - (L - end_decay)) / end_decay * 0.85), 1.0)
    return amp

def signals_from_alleles(alleles, amp, noise_frac=0.01, seed=0, artifact=None, context_ref=None):
    """alleles: list of (sequence, fraction); all sequences equal length L. Returns (L,4) heights and phred-ish scores."""
    rng = np.random.default_rng(seed + 7)
    L = len(alleles[0][0])
    comp = np.zeros((L, 4))
    for seq, frac in alleles:
        assert len(seq) == L
        for j, b in enumerate(seq):
            if b in BASES: comp[j, BASES.index(b)] += frac
            else: comp[j] += frac / 4.0
    if artifact is not None:
        # context-specific minor peaks (shoulders): present wherever the allele's local sequence equals the reference's
        comp = np.zeros((L, 4))
        for seq, frac in alleles:
            for j, b in enumerate(seq):
                row = np.zeros(4)
                if b in BASES: row[BASES.index(b)] = 1.0
                else: row[:] = 0.25
                if j in artifact and j >= 3 and seq[j-3:j+1] == context_ref[j-3:j+1]:
                    sec, a = artifact[j]; row = (1 - a) * row; row[sec] += a
                comp[j] += frac * row
    comp += np.abs(rng.normal(0, noise_frac, comp.shape))
    comp /= comp.sum(axis=1, keepdims=True)
    scale = np.array([CH_SCALE[b] for b in BASES])
    return comp * amp[:, None] * scale[None, :], comp

def write_ab1(path, seq_calls, signals, seed=0, noise=0.012):
    L = signals.shape[0]
    chans, ploc = make_trace(signals, seed=seed, noise=noise)
    order = [BASES.index(b) for b in FWO]
    data = [np.clip(chans[o], 0, 32000).round().astype(">i2").tobytes() for o in order]
    calls = "".join(seq_calls).encode()
    main = signals.max(axis=1) / np.maximum(signals.sum(axis=1), 1)
    qual = np.clip(3 * main * 100 - 240, 0, 60).astype(np.uint8).tobytes()
    ploc_b = ploc.astype(">i2").tobytes()
    ents = [_entry("DATA", 9 + k, 4, 2, len(chans[0]), data[k]) for k in range(4)]
    ents += [_entry("FWO_", 1, 2, 1, 4, FWO.encode()),
             _entry("PBAS", 1, 2, 1, L, calls), _entry("PBAS", 2, 2, 1, L, calls),
             _entry("PCON", 1, 2, 1, L, qual), _entry("PCON", 2, 2, 1, L, qual),
             _entry("PLOC", 1, 4, 2, L, ploc_b), _entry("PLOC", 2, 4, 2, L, ploc_b)]
    open(path, "wb").write(build_abif(sorted(ents, key=lambda e: (e["name"], e["num"]))))
