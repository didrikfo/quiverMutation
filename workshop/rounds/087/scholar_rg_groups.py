"""Round 087 scholar: pick representatives of the periodic-Phi indefinite groups at n (E-190, E-209, E-242 grouping by Coxeter
polynomial) and non-periodic indefinite controls, run the rank-growth test of scholar_rankgrowth.py on each.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/087/scholar_rg_groups.py N K NREP [--ctrl NCTRL] [--cap C] [--only STR,STR]"""
import sys, collections, time
import numpy as np
sys.path.insert(0, "workshop/rounds/087")
import scholar_rankgrowth as rg
from quivermutation import coxeterTables as ct, nakayama as nk
N, K, NREP = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
NC = int(sys.argv[sys.argv.index("--ctrl") + 1]) if "--ctrl" in sys.argv else 0
cap = int(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 300000
only = sys.argv[sys.argv.index("--only") + 1].split(",") if "--only" in sys.argv else None
def order(Phi, cap=400):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all() or (M == -I).all(): return m
        M = M @ Phi
    return 0
groups = collections.defaultdict(list); nonper = []
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    o = order(Phi)
    if o == 0: nonper.append(("".join(map(str, rl)), np.abs(np.linalg.eigvals(Phi)).max()))
    else: groups[tuple(sorted(np.round(np.angle(np.linalg.eigvals(Phi)), 6)))].append(("".join(map(str, rl)), o))
def go(tag, s):
    t = time.time(); A, res, cur = rg.run([int(c) for c in s], K, cap, verbose=False)
    tot = [sum(r) if r else None for r in res]
    print(f"{tag} {s}: totals {tot}  ({time.time()-t:.0f}s)", flush=True)
print(f"n={N}: {len(groups)} periodic groups sizes {sorted(len(v) for v in groups.values())}; {len(nonper)} non-periodic indefinite")
for g, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    mem.sort(key=lambda m: (m[0].count('0'), m[0]))
    reps = [m for m in mem if (only is None or m[0] in only)][:NREP]
    for s, o in reps: go(f"group({len(mem)}) ordPhi={o}", s)
nonper.sort(key=lambda m: m[1])
for s, rho in nonper[:NC]: go(f"CONTROL non-periodic rho(Phi)={rho:.3f}", s)
