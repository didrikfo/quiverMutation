"""Round 071 (idea): finer Z[Phi]-isometry invariant. For each periodic-indefinite n = 11 group member, the canonical sublattices
L_f = ker f(Phi) (f = cyclotomic factor of the Coxeter polynomial; saturated in Z^n) carry the restriction of G = C^-1 + C^-T.
Isometry class of (L_f, G|L_f) is an invariant of the Z[Phi]-isometry class of (Z^n, G, Phi), equivalently of Cartan congruence
(C^-1 = (I-Phi)^-1 G since 1 is not an eigenvalue). Reports det, signature, and the theta series (counts of v with |G(v)| = m,
m <= MAXM) of each L_f; a NO is: two members of one group with different (det, signature, theta) on some L_f.
NULL: identical on all L_f for all members. Control: members of different groups differ (checks the invariant is not blind).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/071/maverick_pieces.py [N] [MAXM]"""
import sys, time, collections, itertools
import numpy as np, sympy
from sympy import Matrix, QQ
from sympy.polys.matrices import DomainMatrix
from math import lcm
sys.path.insert(0, "workshop/rounds/071")
from maverick_conj import refine
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
MAXM = int(sys.argv[2]) if len(sys.argv) > 2 else 8
t0 = time.time()
t = sympy.symbols('t')

def satKernel(M, n):
    """saturated integer basis (rows) of {x in Z^n : M x = 0}."""
    dm = DomainMatrix([[QQ(int(x)) for x in M.row(i)] for i in range(M.rows)], (M.rows, n), QQ)
    R, piv = dm.rref(); R = R.to_Matrix(); free = [c for c in range(n) if c not in piv]
    D = 1
    for i in range(len(piv)):
        for c in free: D = lcm(D, int(sympy.Rational(R[i, c]).q))
    basis = [[1 if i == j else 0 for j in range(len(free))] for i in range(len(free))]
    for i in range(len(piv)):
        row = [int(-R[i, c] * D) for c in free]
        if any(x % D for x in row): basis = refine(basis, row, D)
    out = []
    for w in basis:
        x = [0] * n
        for jj, c in enumerate(free): x[c] = w[jj]
        for i, pc in enumerate(piv): x[pc] = sum(-R[i, c] * w[jj] for jj, c in enumerate(free))
        out.append([int(v) for v in x])
    return out

def gramLLL(B, G, delta=0.75):
    """LLL on rows B w.r.t. positive definite Gram form G (numpy float)."""
    B = [list(r) for r in B]; n = len(B)
    def ip(a, b): return float(np.array(a) @ G @ np.array(b))
    def gs():
        Bs = []; mu = np.zeros((n, n))
        for i in range(n):
            v = np.array(B[i], float)
            for j in range(i):
                mu[i, j] = float(np.array(B[i]) @ G @ Bs[j]) / float(Bs[j] @ G @ Bs[j]); v = v - mu[i, j] * Bs[j]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gs(); k = 1; guard = 0
    while k < n and guard < 20000:
        guard += 1
        for j in range(k - 1, -1, -1):
            q = int(round(mu[k, j]))
            if q: B[k] = [a - q * b for a, b in zip(B[k], B[j])]; Bs, mu = gs()
        if Bs[k] @ G @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ G @ Bs[k - 1]): k += 1
        else: B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(); k = max(k - 1, 1)
    return B

def theta(Gram, maxm, R=4):
    r = Gram.shape[0]; cnt = collections.Counter()
    for c in itertools.product(range(-R, R + 1), repeat=r):
        v = np.array(c); q = int(v @ Gram @ v)
        if 0 < abs(q) <= maxm: cnt[q] += 1
    return tuple(sorted(cnt.items()))

def phiOrd(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return m
        M = M @ Phi
    return 0

keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    if phiOrd(Phi) == 0: continue
    groups[keys[rl]].append((rl, np.rint(Ci + Ci.T).astype(np.int64), Phi))
print(f"groups {sorted(len(v) for v in groups.values())} {time.time()-t0:.0f}s", flush=True)
summary = {}
for key, mem in sorted(groups.items(), key=lambda kv: len(kv[1])):
    rl0, G0, Phi0 = mem[0]
    facs = [f for f, m in sympy.factor_list(Matrix(Phi0.tolist()).charpoly(t).as_expr())[1]]
    print(f"group {len(mem)}: factors {facs}", flush=True)
    invs = collections.defaultdict(collections.Counter)
    for rl, G, Phi in mem:
        PhiM = Matrix(Phi.T.tolist()); row = []
        for f in facs:
            fP = sympy.zeros(N, N); Pw = sympy.eye(N)
            for c in reversed(sympy.Poly(f, t).all_coeffs()):
                fP += int(c) * Pw; Pw = Pw * PhiM
            Bk = satKernel(fP, N)
            Gram = np.array(Bk, dtype=object) @ G.astype(object) @ np.array(Bk, dtype=object).T
            Gf = np.array(Gram, dtype=np.int64); ev = np.linalg.eigvalsh(Gf.astype(float))
            sg = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
            if len(Bk) <= 6 and (sg[0] == 0 or sg[1] == 0):
                sgn = 1 if sg[1] == 0 else -1
                r_ = len(Bk); T = gramLLL([[int(i == j) for j in range(r_)] for i in range(r_)], (sgn * Gf).astype(float))
                Gram2 = np.array(T, dtype=np.int64) @ Gf @ np.array(T, dtype=np.int64).T
                th = theta(sgn * Gram2, MAXM, 3 if len(Bk) > 3 else 5)
            else: th = ("indef/big", round(abs(np.linalg.det(Gf.astype(float)))))
            row.append((len(Bk), sg, round(np.linalg.det(Gf.astype(float))), th))
        invs[tuple(str(r) for r in row)][rl] += 1
    print(f"  distinct (rank, sig, det, theta) tuples over L_f: {len(invs)}", flush=True)
    for k, v in list(invs.items())[:6]:
        print("   ", sum(v.values()), "members:", [x for x in k][:4])
    summary[len(mem)] = set(invs)
print("cross-group disjoint (control):", all(not (summary[a] & summary[b]) for a in summary for b in summary if a < b))
print(f"total {time.time()-t0:.0f}s")
