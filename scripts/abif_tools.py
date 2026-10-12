
import struct, numpy as np
def abif_dir(buf):
    assert buf[:4] == b"ABIF"
    name, num, etype, esize, n, dsize, off, handle = struct.unpack(">4siHHiiii", buf[6:34])
    ents = {}
    for k in range(n):
        o = off + 28*k
        nm, nu, et, es, ne, ds, do, hd = struct.unpack(">4siHHiiii", buf[o:o+28])
        ents[(nm.decode(), nu)] = dict(entry=o, type=et, size=es, num=ne, datasize=ds, offset=(o+20 if ds <= 4 else do))
    return ents
def get_entry(buf, ents, key, dtype=None):
    e = ents[key]; raw = buf[e["offset"]: e["offset"]+e["datasize"]]
    return raw
def patch_abif(src, dst, replacements):
    """replacements: {(tag,num): bytes of identical length}"""
    buf = bytearray(open(src, "rb").read()); ents = abif_dir(bytes(buf))
    for key, new in replacements.items():
        e = ents[key]; assert len(new) == e["datasize"], (key, len(new), e["datasize"])
        buf[e["offset"]: e["offset"]+len(new)] = new
    open(dst, "wb").write(bytes(buf))
def read_channels(path):
    buf = open(path, "rb").read(); ents = abif_dir(buf)
    fwo = get_entry(buf, ents, ("FWO_", 1)).decode()
    ch = {}
    for n, base in zip((9,10,11,12), fwo):
        raw = get_entry(buf, ents, ("DATA", n)); ch[base] = np.frombuffer(raw, dtype=">i2").astype(float)
    ploc = np.frombuffer(get_entry(buf, ents, ("PLOC", 1)), dtype=">i2").astype(int)
    pbas = get_entry(buf, ents, ("PBAS", 1)).decode()
    return dict(fwo=fwo, ch=ch, ploc=ploc, pbas=pbas, ents=ents)
