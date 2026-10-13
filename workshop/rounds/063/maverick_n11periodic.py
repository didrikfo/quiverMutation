"""Round 063 response, item 3: n = 11 check of Candidate C. Among all LNAs at length N (default 11) with periodic Phi and
indefinite Euler form (F-048 rows), group by Coxeter key (all LNAs sharing the key are listed), and report for each group:
size, ord Phi, smallest a with Phi^a = -I, distinct Euler signatures, and the number of members per signature (a group
with >= 2 members of one signature is where Cartan-level C could be live; a group with one member per signature is
already separated by signature). usage (repo root): timeout 10m .venv/bin/python workshop/rounds/063/maverick_n11periodic.py [N] [LIMIT]"""
import sys, time, collections
import numpy as np
from quivermutation import coxeterTables as ct, nakayama as nk

N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
t0 = time.time()

def phiInfo(C, cap=200):
    Ci = np.linalg.inv(C)
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    I = np.eye(len(C), dtype=np.int64); M = Phi.copy(); neg = None
    for m in range(1, cap + 1):
        if neg is None and (M == -I).all(): neg = m
        if (M == I).all(): return m, neg
        M = M @ Phi
        if np.abs(M).max() > 10**9: return 0, None
    return 0, None

def sig(C):
    ev = np.linalg.eigvalsh(C + C.T)
    return (int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()), int((ev > 1e-9).sum()))

allL = list(nk.allRelationLengths(N))[:LIMIT]
keys = ct.lnaKeys(N)
rows = {}
for rl in allL:
    C = np.array(ct.lnaCartanMatrix(N, rl), float)
    Ci = np.linalg.inv(C); G = Ci + Ci.T
    indef = np.linalg.eigvalsh(G).min() < -1e-9
    if not indef: continue
    o, neg = phiInfo(C)
    if o == 0: continue
    rows[rl] = (o, neg, sig(C), tuple(np.rint(np.linalg.eigvalsh(G) * 1).round(6)))
print(f"n={N}: LNAs scanned {len(allL)}; F-048 rows (periodic, indefinite) {len(rows)}; {time.time()-t0:.1f}s")
groups = collections.defaultdict(list)
for rl in rows: groups[keys[rl]].append(rl)
print(f"Coxeter keys among F-048 rows: {len(groups)}")
live = []; sizes = collections.Counter()
for key, mem in groups.items():
    # all LNAs with this key (members need not be F-048 rows if not periodic/indefinite -- same key, so report)
    sizes[len(mem)] += 1
    sigc = collections.Counter(rows[r][2] for r in mem)
    if len(mem) >= 2 and max(sigc.values()) >= 2: live.append((key, mem))
print("group sizes (members: number of keys):", dict(sorted(sizes.items())))
print("ord Phi / smallest a with Phi^a=-I:", dict(collections.Counter((rows[r][0], rows[r][1]) for r in rows)))
print("signatures:", dict(collections.Counter(rows[r][2] for r in rows)))
mixed = [k for k, m in groups.items() if len({rows[r][2] for r in m}) > 1]
print(f"keys with >1 signature among F-048 rows (separated by Euler signature): {len(mixed)}")
print(f"keys with >=2 F-048 members of one signature (Cartan-level C potentially live): {len(live)}")
for key, mem in live[:60]:
    print("  ", [ct.className(r) for r in mem], rows[mem[0]][:3])
# also: do any F-048 rows share a key with an LNA that is NOT an F-048 row?
other = sum(1 for k in groups for r in ct.lnaKeyIndex(N)[k] if r not in rows)
print("LNAs sharing a key with an F-048 row but not themselves F-048 rows:", other)
print(f"total {time.time()-t0:.1f}s")
