"""Round 074 (idea): power survey of maverick_isometry.isometry at n = 10 over ALL definite-G LNAs with periodic Phi (cap ord 200).
Group by Coxeter key; inside a key group test each member against the first (YES/NO/other); across groups with equal charpoly+detG test reps.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/074/maverick_survey10.py [N=10]"""
import sys, time, collections
sys.path.insert(0, "workshop/rounds/074")
import numpy as np, sympy
from maverick_isometry import data, isometry
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10; t0 = time.time()
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
def ordOf(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return m
        M = M @ Phi
    return 0
tot = 0
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, tuple(rl)), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() < 1e-9: continue
    tot += 1
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    if ordOf(Phi) == 0: continue
    groups[keys[rl]].append(rl)
print(f"definite LNAs {tot}; periodic definite groups {len(groups)} sizes {sorted(len(v) for v in groups.values())} {time.time()-t0:.0f}s", flush=True)
res = collections.Counter(); reps = []
for k, mem in groups.items():
    a = mem[0]; Ga, Ma, _ = data(a, N); reps.append((k, a, Ga, Ma))
    c = collections.Counter()
    for b in mem[1:]:
        Gb, Mb, _ = data(b, N)
        r, _, _ = isometry(Ga, Ma, Gb, Mb); c[r] += 1
        if time.time() - t0 > 480: break
    print("group", len(mem), "rep", "".join(map(str, a)), dict(c), f"{time.time()-t0:.0f}s", flush=True); res.update(c)
print("within-key totals", dict(res))
cross = collections.Counter()
for i in range(len(reps)):
    for j in range(i + 1, len(reps)):
        ki, a, Ga, Ma = reps[i]; kj, b, Gb, Mb = reps[j]
        if round(np.linalg.det(Ga)) != round(np.linalg.det(Gb)): cross["detG differs"] += 1; continue
        r, _, _ = isometry(Ga, Ma, Gb, Mb); cross[r] += 1
        if r != "NO": print("cross", "".join(map(str, a)), "".join(map(str, b)), r)
print("cross-key rep pairs", dict(cross), f"{time.time()-t0:.0f}s")
