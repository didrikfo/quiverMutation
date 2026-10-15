"""Round 066: Candidate C at n = 11. For the periodic-indefinite LNA groups (E-190), compute Cartan-level finer invariants
beyond signature: Smith form of the integral symmetrised Euler form G = C^-1 + C^-T (congruence invariant under
derived equivalence, since C' = P^T C P with P in GL_n(Z)), Smith form of C (trivial, det 1), and of Phi - I, Phi + I,
Phi^2 + Phi + 1 etc. (Coxeter-polynomial factor lattices). Counts distinct values per group.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/066/maverick_finer.py [N] [LIMIT]"""
import sys, time, collections
import numpy as np, sympy
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
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
allL = list(nk.allRelationLengths(N))[:LIMIT]
for rl in allL:
    C = np.array(ct.lnaCartanMatrix(N, rl), float)
    Ci = np.linalg.inv(C); G = Ci + Ci.T
    if np.linalg.eigvalsh(G).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    if phiOrd(Phi) == 0: continue
    groups[keys[rl]].append((rl, np.rint(G).astype(np.int64), Phi))
print(f"n={N}: groups {len(groups)}, sizes {sorted(len(v) for v in groups.values())} {time.time()-t0:.1f}s")
for key, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    c = collections.Counter(); cphi = collections.Counter(); cdet = collections.Counter()
    ex = {}
    I = np.eye(N, dtype=np.int64)
    for rl, G, Phi in mem:
        s = snf(G)
        s2 = (s, snf(Phi - I), snf(Phi + I))
        c[s] += 1; cphi[s2] += 1; cdet[round(np.linalg.det(G))] += 1
        ex.setdefault(s2, ct.className(rl))
    print(f"group size {len(mem)}: ord Phi {phiOrd(mem[0][2])}; det G {dict(cdet)}; SNF(G) classes {len(c)} {dict(c)}; "
          f"SNF(G,Phi-I,Phi+I) classes {len(cphi)} {sorted(cphi.values(), reverse=True)}")
    for s2, e in list(ex.items())[:3]: print("    ex", e, s2)
print(f"total {time.time()-t0:.1f}s")
