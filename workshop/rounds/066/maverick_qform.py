"""Round 066: discriminant quadratic form q_G of the Euler lattice (Z^n, G), G = C^-1 + C^-T (even, diag 2), on A = coker G.
For x = G^-1 y in L*/L, q(x) = y^T G^-1 y mod 2. Invariant per member: multiset of (order of x, q(x)) over ALL x in A
(A has order det G <= 8); congruence-invariant. Groups = same as maverick_finer.py (indefinite, periodic Phi, same key).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/066/maverick_qform.py [N]"""
import sys, time, collections
from fractions import Fraction as Fr
import numpy as np
from sympy import Matrix, ZZ, QQ
from sympy.polys.matrices import DomainMatrix
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
t0 = time.time()
def phiPeriodic(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return True
        M = M @ Phi
    return False
def qinv(G):
    Gi = DomainMatrix.from_Matrix(Matrix(G.tolist())).convert_to(QQ).inv().to_Matrix()
    Gi = [[Fr(int(Gi[i, j].p), int(Gi[i, j].q)) for j in range(N)] for i in range(N)]
    def red(v): return tuple(x - (x.numerator // x.denominator) for x in v)   # mod Z^n
    start = [red(tuple(Gi[i][j] for i in range(N))) for j in range(N)]       # classes G^-1 e_j
    zero = tuple(Fr(0) for _ in range(N)); elems = {zero}; frontier = [zero]
    while frontier:
        nf = []
        for e in frontier:
            for s in start:
                t = red(tuple(a + b for a, b in zip(e, s)))
                if t not in elems: elems.add(t); nf.append(t)
        frontier = nf
    out = collections.Counter()
    for x in elems:
        # x = G^-1 y with y = G x (integer); q = x^T G x mod 2
        Gx = [sum(int(G[i, j]) * x[j] for j in range(N)) for i in range(N)]
        q = sum(x[i] * Gx[i] for i in range(N)) % 2
        o = 1
        while any((o * a).denominator != 1 for a in x): o += 1
        out[(o, str(q))] += 1
    return tuple(sorted(out.items())), len(elems)
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C); G = Ci + Ci.T
    if np.linalg.eigvalsh(G).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    if not phiPeriodic(Phi): continue
    groups[keys[rl]].append(np.rint(G).astype(np.int64))
print(f"n={N}: groups {len(groups)} sizes {sorted(len(v) for v in groups.values())} {time.time()-t0:.1f}s")
for k, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    c = collections.Counter(); 
    for G in mem:
        assert (np.diag(G) == 2).all()
        r, a = qinv(G); c[r] += 1
    print(f"group size {len(mem)} |A|={a} det {round(np.linalg.det(mem[0]))}: distinct q-forms {len(c)} counts {sorted(c.values(), reverse=True)}")
    for r in c: print("   (order,q mod 2):count", r)
print(f"total {time.time()-t0:.1f}s")
