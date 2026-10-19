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

def orbitCols(M, vs, n):
    """greedily pick columns M^a v_i (a ascending, round robin) until rank n; return list of (i, a) and the matrix (object)."""
    cols = []; idx = []; rk = 0
    cur = [v.astype(object) for v in vs]
    for a in range(n):
        for i in range(len(vs)):
            trial = cols + [cur[i]]
            if np.linalg.matrix_rank(np.array(trial, dtype=float)) > rk:
                cols = trial; idx.append((i, a)); rk += 1
        cur = [M.astype(object) @ c for c in cur]
        if rk == n: break
    return idx, (Matrix(np.array(cols, dtype=object).T.tolist()) if rk == n else None)

def cols_for(M, us, idx):
    out = []
    for (i, a) in idx:
        w = us[i].astype(object)
        for _ in range(a): w = M.astype(object) @ w
        out.append(w)
    return Matrix(np.array(out, dtype=object).T.tolist())

def isometry(Ga, Ma, Gb, Mb, cap=300000, kmax=4):
    """Pick generators v_1..v_k of least norm until the M-orbits span Q^n; enumerate images u_i level by level (norm AND all
    G(v_i, M^c v_j) = G(u_i, M^c u_j) pruned); P = W(u) W(v)^-1. Exhaustive => NO is a proof of non-isometry (within GL_n(Z) equivariant)."""
    n = Ga.shape[0]
    if np.linalg.eigvalsh(Ga.astype(float)).min() < 1e-9: return "SKIP-indef", None, 0
    Mp = [np.linalg.matrix_power(Ma, c) for c in range(n)]; Mq = [np.linalg.matrix_power(Mb, c) for c in range(n)]
    def span(ws): return np.linalg.matrix_rank(np.array([Mp[c] @ w for w in ws for c in range(n)], dtype=float)) if ws else 0
    vs = []
    for v in sorted([v for v in shortVectors(Ga, 16 if n <= 4 else 4) if v.any()], key=lambda z: int(z @ Ga @ z)):
        if span(vs + [v]) > span(vs): vs = vs + [v]
        if span(vs) == n: break
    if span(vs) < n or len(vs) > kmax: return "SKIP-gens", None, 0
    idx, W = orbitCols(Ma, vs, n)
    Winv = W.inv()
    Bmax = max(int(v @ Ga @ v) for v in vs)
    byNorm = collections.defaultdict(list)
    for u in shortVectors(Gb, Bmax):
        if u.any(): byNorm[int(u @ Gb @ u)].append(u)
    nodes = [0]; found = [None]
    def rec(k, us):
        if found[0] is not None or nodes[0] > cap: return
        if k == len(vs):
            P = cols_for(Mb, us, idx) * Winv
            if any(x.q != 1 for x in P) or abs(P.det()) != 1: return
            Ga_m = Matrix(Ga.tolist()); Gb_m = Matrix(Gb.tolist())
            if P.T * Gb_m * P == Ga_m and P * Matrix(Ma.tolist()) == Matrix(Mb.tolist()) * P: found[0] = P
            return
        q = int(vs[k] @ Ga @ vs[k])
        for u in byNorm[q]:
            nodes[0] += 1
            ok = all(int(vs[k] @ Ga @ Mp[c] @ vs[j]) == int(u @ Gb @ Mq[c] @ us[j]) and int(vs[j] @ Ga @ Mp[c] @ vs[k]) == int(us[j] @ Gb @ Mq[c] @ u)
                     for j in range(k) for c in range(n)) and all(int(vs[k] @ Ga @ Mp[c] @ vs[k]) == int(u @ Gb @ Mq[c] @ u) for c in range(n))
            if ok: rec(k + 1, us + [u])
    rec(0, [])
    if found[0] is not None: return "YES", found[0], nodes[0]
    return ("INCONCLUSIVE-cap" if nodes[0] > cap else "NO"), None, nodes[0]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    t0 = time.time()
    rd = lambda s: tuple(int(c) for c in s)
    A, B = rd("90000000"), rd("00000030")
    for a, b in [(A, B), (B, A), (A, A)]:
        Ga, Ma, _ = data(a, N); Gb, Mb, _ = data(b, N)
        res, P, k = isometry(Ga, Ma, Gb, Mb)
        print(a, b, res, "candidates", k, f"{time.time()-t0:.1f}s", flush=True)
