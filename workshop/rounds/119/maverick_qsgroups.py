"""Round 119 maverick: do the nonsingular n = 12 low-lam lnaKeys groups (lam < 1.5; |det S| = 3, 4, 8, 11) contain members with DIFFERENT discriminant forms q_S?
q_S fingerprint = sorted histogram of x^T S x mod 2 over x in S^-1 Z^n / Z^n (fractions, exact).  Different fingerprint => no unimodular P with P C_a P^T = C_b (Nikulin 1.9.4 on S = C + C^T).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_qsgroups.py [members per group = 400]"""
import sys, collections, random; sys.path.insert(0, "workshop/rounds/097")
import numpy as np, sympy as sp
from fractions import Fraction
from maverick_cartanZ import cartan, disc_form
from quivermutation import coxeterTables as ct
N = 12; per = int(sys.argv[1]) if len(sys.argv) > 1 else 400; random.seed(1)
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, kk in keys.items(): groups[kk].append(tuple(rl))
def lam_of(rl):
    C = np.array(ct.lnaCartanMatrix(N, list(rl))).astype(float)
    if round(np.linalg.det(C)) != 1: return None
    ev = np.linalg.eigvalsh(C + C.T)
    if (ev < -1e-9).sum() != 1 or (abs(ev) < 1e-9).any(): return None
    w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C); big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-6]
    return big[0] if len(big) == 1 else None
def fp(rl):
    S = cartan(N, list(rl)); S = S + S.T; d, E, G, bil, q, add = disc_form(S)
    Sf = [[Fraction(int(S[i, j])) for j in range(N)] for i in range(N)]
    h = collections.Counter(sum(x[i] * Sf[i][j] * x[j] for i in range(N) for j in range(N)) % 2 for x in E)
    return tuple(sorted(h.items()))
sel = sorted([(lam_of(m[0]), kk, m) for kk, m in groups.items() if lam_of(m[0]) and lam_of(m[0]) < 1.5], key=lambda t: t[0])
for gi, (lam, kk, mem) in enumerate(sel):
    s = random.sample(mem, min(per, len(mem))); hist = collections.defaultdict(list)
    for rl in s: hist[fp(rl)].append(rl)
    print("GROUP", gi, "lam %.6f" % lam, "size", len(mem), "sampled", len(s), "fingerprints", {str(k): len(v) for k, v in hist.items()}, flush=True)
    if len(hist) > 1:
        for k, v in hist.items(): print("   ", k, v[:3], flush=True)
