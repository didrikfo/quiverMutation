"""Round 074 (idea): staged explicit-isometry test (Cartan congruence = Phi^T-equivariant isometry of G; E-209 convention).
Stage 1 (this script, exact, definite G only): for (G_a, M_a=Phi_a^T), (G_b, M_b): pick v with W(v)=[v, Mv, ..] invertible over Q and
minimal G_a(v); every equivariant isometry P has P v = u with G_b(u) = G_a(v), P = W(u) W(v)^-1. Enumerate ALL such u (Fincke-Pohst), test
P integral, det +-1, P^T G_b P = G_a, P M_a = M_b P. 'YES' = explicit P (exact int check). 'NO' = exhaustive failure over every candidate u
(sound only if Q[M_a] has cyclic vector v, checked). 'SKIP' = G indefinite or no cyclic vector found.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/074/maverick_isometry.py [N=10] """
import sys, time, itertools, collections
import numpy as np
from fractions import Fraction
import sympy
from sympy import Matrix
from quivermutation import coxeterTables as ct, nakayama as nk

def data(rl, N):
    C = np.array(ct.lnaCartanMatrix(N, tuple(rl)), float); Ci = np.linalg.inv(C)
    Phi = np.rint(-Ci.T @ C).astype(np.int64); G = np.rint(Ci + Ci.T).astype(np.int64)
    M = Phi.T.copy()
    assert (M.T @ G @ M == G).all(), "M does not preserve G"
    return G, M, Phi

def shortVectors(H, B):
    """all x in Z^r with x^T H x <= B (H pos. def. integer). Fincke-Pohst on the Cholesky of H (no LLL: fine for r <= 11, small B)."""
    r = H.shape[0]; Hf = H.astype(float); L = np.linalg.cholesky(Hf).T   # H = L^T L, L upper
    out = []; x = [0] * r
    def rec(i, rem):
        if i < 0:
            out.append(tuple(x)); return
        s = sum(L[i, j] * x[j] for j in range(i + 1, r))
        rad = max(rem, 0) ** 0.5 / L[i, i]; c = -s / L[i, i]
        for xi in range(int(np.ceil(c - rad - 1e-9)), int(np.floor(c + rad + 1e-9)) + 1):
            x[i] = xi; t = (L[i, i] * xi + s) ** 2
            if t <= rem + 1e-7: rec(i - 1, rem - t)
        x[i] = 0
    rec(r - 1, B + 1e-7)
    return [np.array(v, dtype=np.int64) for v in out]

def Wmat(M, v, n):
    cols = [v.astype(object)]
    for _ in range(n - 1): cols.append(M.astype(object) @ cols[-1])
    return Matrix(np.array(cols, dtype=object).T.tolist())

def isometry(Ga, Ma, Gb, Mb, cap=2_000_000):
    n = Ga.shape[0]
    if np.linalg.eigvalsh(Ga.astype(float)).min() < 1e-9: return "SKIP-indef", None, 0
    if sorted(np.linalg.eigvalsh(Ga.astype(float)).round(6)) != sorted(np.linalg.eigvalsh(Gb.astype(float)).round(6)) and False: pass
    # cyclic v of least norm: scan short vectors of Ga by norm
    best = None
    for B in (2, 4, 6, 8, 12):
        for v in sorted(shortVectors(Ga, B), key=lambda z: int(z @ Ga @ z)):
            if not v.any(): continue
            Wf = np.array([np.linalg.matrix_power(Ma.astype(float), k) @ v for k in range(n)]).T
            if abs(np.linalg.det(Wf)) < 1e-6: continue
            W = Wmat(Ma, v, n)
            if W.det() != 0: best = (v, W); break
        if best: break
    if best is None: return "SKIP-nocyc", None, 0
    v, W = best; q = int(v @ Ga @ v); Winv = W.inv()
    cands = [u for u in shortVectors(Gb, q) if int(u @ Gb @ u) == q]
    if len(cands) > cap: return "SKIP-big", None, len(cands)
    Ga_m = Matrix(Ga.tolist()); Gb_m = Matrix(Gb.tolist()); Ma_m = Matrix(Ma.tolist()); Mb_m = Matrix(Mb.tolist())
    for u in cands:
        P = Wmat(Mb, u, n) * Winv
        if any(x.q != 1 for x in P): continue
        if abs(P.det()) != 1: continue
        if P.T * Gb_m * P == Ga_m and P * Ma_m == Mb_m * P: return "YES", P, len(cands)
    return "NO", None, len(cands)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    t0 = time.time()
    rd = lambda s: tuple(int(c) for c in s)
    A, B = rd("90000000"), rd("00000030")
    for a, b in [(A, B), (B, A), (A, A)]:
        Ga, Ma, _ = data(a, N); Gb, Mb, _ = data(b, N)
        res, P, k = isometry(Ga, Ma, Gb, Mb)
        print(a, b, res, "candidates", k, f"{time.time()-t0:.1f}s", flush=True)
