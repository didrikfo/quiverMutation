"""Round 071 (idea): Z-conjugacy of Phi inside the n = 11 groups (F-047 level, i.e. Z[Phi]-module iso of (Z^n, Phi)).
For members a, b of one group: find P in GL_n(Z) with P Phi_a = Phi_b P. Method: Phi_a has a cyclic-ish vector v (Z{Phi^k v}
of finite index d); P is fixed by w = P v; integrality of P = W(w) V^-1 is a congruence system on w mod d, solved by lattice
refinement; then LLL-reduce the lattice Lambda of admissible w and search small combinations for det P = +-1.
YES = explicit P found (exact integer check). NO would be: Phi_a, Phi_b separated by an invariant (see maverick_snfx.py)
or an exhaustive failure -- this script cannot prove NO, so 'not found' is INCONCLUSIVE, reported as such.
Controls: (c1) b = Q^-1 Phi_a Q for random unimodular Q (must be found); (c2) Phi restricted to an index-3/other
Phi-invariant sublattice (not Z-conjugate if SNF probes differ; script must NOT find P).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/071/maverick_conj.py [N] [GROUP_SIZE] [TRIALS]"""
import sys, time, random, collections, itertools
from fractions import Fraction
import numpy as np, sympy
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct, nakayama as nk

def phiOf(rl, N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    return np.rint(-Ci.T @ C).astype(np.int64).tolist()

def matmul(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def lll(B, delta=0.75):
    """float LLL on rows of an integer matrix (list of lists); integer rows out (unimodular transform applied exactly)."""
    B = [list(map(int, r)) for r in B]; n = len(B)
    def gs():
        A = np.array(B, dtype=float); Bs = np.zeros_like(A); mu = np.zeros((n, n))
        for i in range(n):
            v = A[i].copy()
            for j in range(i):
                mu[i, j] = A[i] @ Bs[j] / (Bs[j] @ Bs[j]); v -= mu[i, j] * Bs[j]
            Bs[i] = v
        return Bs, mu
    Bs, mu = gs(); k = 1; guard = 0
    while k < n and guard < 100000:
        guard += 1
        for j in range(k - 1, -1, -1):
            q = int(round(mu[k, j]))
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]; Bs, mu = gs()
        if Bs[k] @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ Bs[k - 1]): k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(); k = max(k - 1, 1)
    return B

def egcd(a, b):
    if b == 0: return (abs(a), (1 if a >= 0 else -1), 0)
    g, x, y = egcd(b, a % b); return g, y, x - (a // b) * y

def refine(basis, row, d):
    """sublattice of span(basis) (rows) where sum_j row[j]*w_j == 0 mod d."""
    k = len(basis); a = [sum(row[j] * basis[i][j] for j in range(len(row))) % d for i in range(k)]
    # kernel of Z^k -> Z/d, c -> sum c_i a_i. Use gcd chain.
    new = []; g = d; coef = None
    # find unimodular transform: process pairs
    vecs = [list(b) for b in basis]; vals = a[:]
    # reduce to one generator with value g0 = gcd(vals, d), others value 0
    idx = [i for i in range(k) if vals[i] % d]
    if not idx: return vecs
    p = idx[0]
    for i in range(k):
        if i == p or vals[i] % d == 0: continue
        g, x, y = egcd(vals[p], vals[i])
        # new p = x*p + y*i (value g), new i = (vals[i]/g)*p - (vals[p]/g)*i (value 0)
        vp, vi = vecs[p], vecs[i]
        np_ = [x * s + y * t for s, t in zip(vp, vi)]
        ni = [(vals[i] // g) * s - (vals[p] // g) * t for s, t in zip(vp, vi)]
        vecs[p], vecs[i] = np_, ni; vals[p], vals[i] = g, 0
    # p has value g; need multiple m*p with m*g == 0 mod d, m = d/gcd(g,d)
    g0 = egcd(vals[p], d)[0]; m = d // g0
    vecs[p] = [m * x for x in vecs[p]]
    return vecs

def homBasis(PhiA, PhiB, n):
    """Z-basis (list of integer n x n lists) of {P : P PhiA = PhiB P}, via Q-rref then congruence refinement."""
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    rows = []
    for i in range(n):
        for j in range(n):
            r = [0] * (n * n)
            for k in range(n): r[i * n + k] += PhiA[k][j]      # (P PhiA)_ij = sum_k P_ik PhiA_kj
            for k in range(n): r[k * n + j] -= PhiB[i][k]      # (PhiB P)_ij = sum_k PhiB_ik P_kj
            rows.append(r)
    dm = DomainMatrix([[QQ(int(x)) for x in r] for r in rows], (n * n, n * n), QQ)
    R, piv = dm.rref()
    R = R.to_Matrix(); free = [c for c in range(n * n) if c not in piv]
    f = len(free)
    # solution: x_free = e_j ; x_piv_i = -R[i, free_j]
    from math import lcm
    D = 1
    for i in range(len(piv)):
        for jj, c in enumerate(free): D = lcm(D, int(sympy.Rational(R[i, c]).q))
    basis = [[1 if i == j else 0 for j in range(f)] for i in range(f)]
    for i in range(len(piv)):
        row = [int(-R[i, c] * D) for c in free]
        if any(x % D for x in row): basis = refine(basis, row, D)
    out = []
    for w in basis:
        x = [0] * (n * n)
        for jj, c in enumerate(free): x[c] = w[jj]
        for i, pc in enumerate(piv): x[pc] = sum(-R[i, c] * w[jj] for jj, c in enumerate(free))
        out.append([[int(x[i * n + j]) for j in range(n)] for i in range(n)])
    return out

def conjugators(PhiA, PhiB, N, rng, tries=20000):
    n = N; Hs = homBasis(PhiA, PhiB, n); r = len(Hs)
    flat = lll([[x for row in H for x in row] for H in Hs])
    Pf = [np.array(f, dtype=float).reshape(n, n) for f in flat]; Pi = [Matrix(np.array(f, dtype=object).reshape(n, n).tolist()) for f in flat]
    A = Matrix(PhiA); Bm = Matrix(PhiB); found = None; t = 0
    for t in range(tries):
        c = [rng.choice([-1, 0, 0, 1]) if t > 2 * r else (1 if i == t % r else 0) for i in range(r)]
        if not any(c): continue
        P = sum(ci * Pk for ci, Pk in zip(c, Pf))
        if abs(abs(np.linalg.det(P)) - 1) < 1e-6:
            Pm = sum((ci * Pk for ci, Pk in zip(c, Pi)), sympy.zeros(n, n))
            if abs(Pm.det()) == 1 and Pm * A == Bm * Pm: found = Pm; break
    return found, r, t

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    GS = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    TR = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
    rng = random.Random(1); t0 = time.time()
    keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
    for rl in nk.allRelationLengths(N):
        C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
        if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
        groups[keys[rl]].append(rl)
    grp = [m for m in groups.values() if len(m) == GS][0]
    Phis = [phiOf(rl, N) for rl in grp]
    print("group", GS, [ct.className(r) for r in grp][:5], "...", flush=True)
    # control 1: random unimodular conjugate
    A = Phis[0]
    Q = Matrix.eye(N)
    for _ in range(30):
        i, j = rng.sample(range(N), 2); E = Matrix.eye(N); E[i, j] = rng.choice([-1, 1]); Q = Q * E
    B1 = [[int(x) for x in r] for r in (Q.inv() * Matrix(A) * Q).tolist()]
    f, d, t = conjugators(A, B1, N, rng, TR)
    print(f"control1 (random conjugate): {'FOUND' if f is not None else 'not found'} d={d} tries={t} {time.time()-t0:.0f}s", flush=True)
    for b in range(1, len(grp)):
        f, d, t = conjugators(A, Phis[b], N, rng, TR)
        print(f"{ct.className(grp[0])} ~ {ct.className(grp[b])}: {'FOUND' if f is not None else 'not found'} d={d} tries={t} {time.time()-t0:.0f}s", flush=True)
