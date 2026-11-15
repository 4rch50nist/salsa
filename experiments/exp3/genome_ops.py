import sys
from collections import Counter

IDS, BARCODES, PAF, BED, OUT = sys.argv[1:6]
SPAN = 1 << 15

order = [l.strip() for l in open(IDS) if l.strip()]

bc = {}
with open(BARCODES) as f:
    head = f.readline().rstrip("\n").split("\t")
    ri, bi = head.index("read_id"), head.index("barcode")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[bi] != "unclassified":
            bc[p[ri]] = p[bi]

aln = {}
with open(PAF) as f:
    for line in f:
        p = line.split("\t")
        if "tp:A:P" not in p[12:]:
            continue
        rid = p[0]
        if rid not in aln:
            aln[rid] = (int(p[7]), int(p[8]))

left, right = {}, {}
for line in open(BED):
    p = line.split("\t")
    n, side = int(p[3].split("_")[1]), p[3].split("_")[2]
    s, e = int(p[1]), int(p[2])
    if side == "LEFT":
        left[n] = max(left.get(n, 0), e)
    else:
        right[n] = min(right.get(n, 1 << 30), s)
amps = sorted(left)
insert = {n: (left[n], right[n]) for n in amps}
unique = {}
for i, n in enumerate(amps):
    s, e = insert[n]
    if i > 0: s = max(s, insert[amps[i - 1]][1])
    if i + 1 < len(amps): e = min(e, insert[amps[i + 1]][0])
    unique[n] = (s, e)

def amplicon(ts, te):
    best, ov = None, 0
    for n in amps:
        s, e = insert[n]
        o = min(te, e) - max(ts, s)
        if o > ov: best, ov = n, o
    return best

kept = [r for r in order if r in bc and r in aln]
rank = {b: i for i, (b, _) in enumerate(Counter(bc[r] for r in kept).most_common())}

ops = bytearray()
used = 0
for r in kept:
    ts, te = aln[r]
    n = amplicon(ts, te)
    if n is None or te + 1 >= SPAN:
        continue
    s = rank[bc[r]]
    us, ue = unique[n]
    for kind, key in ((2, ts + 1), (2, te + 1), (0, ue), (0, us)):
        if key == 0:
            continue
        ops += bytes((kind, s)) + key.to_bytes(2, "little")
    used += 1
open(OUT, "wb").write(ops)

counts = Counter(bc[r] for r in kept)
top = counts.most_common(1)[0][0]
open("genome_meta.txt", "w").write(f"{len(rank)}\n{top.split('_')[-1]}\n")
print(f"reads selected {len(order):,}, barcoded {sum(r in bc for r in order):,}, aligned {sum(r in aln for r in order):,}, replayed {used:,}")
print(f"samples {len(rank)}, busiest {counts.most_common(1)[0][0]} ({counts.most_common(1)[0][1]:,} reads), operations {len(ops) // 4:,}")
