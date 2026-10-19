"""Round 074 (idea): stage 1 on the n = 11 15-group (E-190, E-209): exact Phi^T-equivariant isometry of the DEFINITE canonical pieces
L_f = ker f(Phi^T) (f = t+1 rank 3, t^2-t+1 rank 2) of all 15 members against member 0, and the 9 orbit+mirror reps. The rank-6 piece
(t^6-...+1, signature (4,2)) is indefinite: SKIPPED (stage 2, not run). YES/NO per piece; a NO would separate the two members at the
Z[Phi]-isometry level (necessary condition for Cartan congruence). Gluing between pieces is NOT tested.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/074/maverick_pieces15.py"""
import sys, time, collections
sys.path.insert(0, "workshop/rounds/074"); sys.path.insert(0, "workshop/rounds/071")
import numpy as np, sympy
from sympy import Matrix, QQ
from sympy.polys.matrices import DomainMatrix
from math import lcm
from maverick_conj import refine
from maverick_isometry import data, isometry
from quivermutation import coxeterTables as ct, nakayama as nk, freeMoves as fm
N = 11; t0 = time.time(); t = sympy.symbols('t')
def satKernel(M, n):
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
def periodic(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return True
        M = M @ Phi
    return False
lnas, orbits = fm.derivedOrbits(N, free=True); parent = {}
for k, mem in orbits.items():
    for m in mem: parent[tuple(m)] = tuple(mem[0])
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for m in list(parent):
    a, b = find(m), find(fm.mirrorRow(N, m))
    if a != b: parent[a] = b
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl in lnas:
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    if not periodic(np.rint(-Ci.T @ C).astype(np.int64)): continue
    groups[keys[rl]].append(tuple(rl))
mem = next(v for v in groups.values() if len(v) == 15)
print("15-group", len(mem), "classes", len({find(m) for m in mem}), f"{time.time()-t0:.0f}s", flush=True)
facs = [t + 1, t**2 - t + 1]
def piece(rl, f):
    G, M, _ = data(rl, N); n = N
    fP = sympy.zeros(n, n); Pw = sympy.eye(n); Mm = Matrix(M.tolist())
    for c in reversed(sympy.Poly(f, t).all_coeffs()):
        fP += int(c) * Pw; Pw = Pw * Mm
    Bk = np.array(satKernel(fP, n), dtype=object).T            # columns = basis
    Gf = np.array(Matrix(Bk.tolist()).T * Matrix(G.tolist()) * Matrix(Bk.tolist()), dtype=np.int64)
    B = Matrix(Bk.tolist()); Mf = (B.T * B).inv() * B.T * Mm * B
    assert all(x.q == 1 for x in Mf)
    return Gf, np.array(Mf.tolist(), dtype=np.int64)
base = [piece(mem[0], f) for f in facs]
tally = collections.Counter(); classOf = {}
for rl in mem:
    row = []
    for f, (Ga, Ma) in zip(facs, base):
        Gb, Mb = piece(rl, f); r, P, nodes = isometry(Ga, Ma, Gb, Mb); row.append(r)
    tally[tuple(row)] += 1
    print("".join(map(str, rl)), "class", "".join(map(str, find(rl))), row, f"{time.time()-t0:.0f}s", flush=True)
print("tally", dict(tally))
