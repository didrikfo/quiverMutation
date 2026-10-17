"""Round 071 (idea): n = 11 periodic-indefinite groups (E-190/E-197). Per group, SNF(x C + C^T) for x in -X..X
(congruence invariant of C -> P^T C P; = SNF(x I - Phi^-1) up to a unit) and SNF(a Phi + b I) for small (a,b);
count distinct values per group. NO-outcome: some group has >= 2 distinct values for some (x) or (a,b) => that group
is NOT one derived class (invariants of Z[Phi]-module + C). NULL: one value per group for every probe.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/071/maverick_snfx.py [N] [X]"""
import sys, time, collections
import numpy as np, sympy
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
X = int(sys.argv[2]) if len(sys.argv) > 2 else 5
t0 = time.time()
def snf(M):
    S = smith_normal_form(sympy.Matrix(M.tolist()), domain=sympy.ZZ)
    return tuple(sorted(abs(int(S[i, i])) for i in range(S.shape[0])))
def phiOrd(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return m
        M = M @ Phi
    return 0
keys = ct.lnaKeys(N)
groups = collections.defaultdict(list)
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float)
    Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    if phiOrd(Phi) == 0: continue
    groups[keys[rl]].append((rl, np.rint(C).astype(np.int64), Phi))
print(f"n={N}: groups sizes {sorted(len(v) for v in groups.values())} {time.time()-t0:.1f}s", flush=True)
I = np.eye(N, dtype=np.int64)
probes = [("x", x) for x in range(-X, X + 1)] + [("ab", ab) for ab in [(2, 1), (3, 1), (3, 2), (2, -1), (3, -1), (4, 1), (5, 2)]]
for key, mem in sorted(groups.items(), key=lambda kv: len(kv[1])):
    res = {}
    for kind, p in probes:
        c = collections.Counter()
        for rl, C, Phi in mem:
            M = (p * C + C.T) if kind == "x" else (p[0] * Phi + p[1] * I)
            c[snf(M)] += 1
        res[(kind, p)] = c
    print(f"group size {len(mem)} ord Phi {phiOrd(mem[0][2])}: distinct SNF counts per probe:")
    print("   ", {f"{k[0]}{k[1]}": len(v) for k, v in res.items()})
    for k, v in res.items():
        if len(v) > 1: print("    SPLIT", k, dict(v))
    print("    x=2:", list(res[("x", 2)].items())[:2], "x=0:", list(res[("x", 0)].items())[:1])
print(f"total {time.time()-t0:.1f}s")
