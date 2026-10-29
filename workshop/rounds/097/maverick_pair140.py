"""Round 097 maverick: the one profile-equal pair NOT separated by q_S (n = 12, key group of 140): Phi data and a module-isomorphism search.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/097/maverick_pair140.py"""
import sys; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
N = 12
a = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); b = (2, 2, 3, 0, 3, 0, 2, 2, 3, 0)
x = sp.symbols('x')
Ca, Cb = cartan(N, a), cartan(N, b)
for C in (Ca, Cb):
    P = phi(C); print(sp.factor(P.charpoly(x).as_expr()), "det C", C.det(), "det S", (C + C.T).det(), "|1-Phi| det", (sp.eye(N) - P).det())
Pa, Pb = phi(Ca), phi(Cb)
print("module-iso random search:", module_iso_search(Pa, Pb, 20000)[:3])

# local test: is there P in the saturated commutant lattice with det P a unit at p?  (exhaustive mod 2, random mod 3,5,7)
import numpy as np, random
n = N
A = np.kron(np.array(Pa.T.tolist(), dtype=object), np.eye(n, dtype=object)) - np.kron(np.eye(n, dtype=object), np.array(Pb.tolist(), dtype=object))
m = n * n
K = pari.matkerint(pari.matrix(m, m, [int(A[i, j]) for i in range(m) for j in range(m)])); k = int(pari.matsize(K)[1])
B = [np.array([int(K[i, j]) for i in range(m)], dtype=np.int64).reshape(n, n) for j in range(k)]
def rank_mod(M, p):
    M = [list(map(int, r)) for r in (M % p)]; r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; iv = pow(M[r][c], -1, p)
        M[r] = [(v * iv) % p for v in M[r]]
        for i in range(n):
            if i != r and M[i][c] % p: f = M[i][c]; M[i] = [(vi - f * vr) % p for vi, vr in zip(M[i], M[r])]
        r += 1
    return r
print("commutant rank", k)
for p, cnt in ((2, None), (3, 20000), (5, 20000), (7, 20000)):
    best = 0; seen = 0
    if cnt is None:
        for bits in range(1 << k):
            M = sum((B[j] for j in range(k) if bits >> j & 1), np.zeros((n, n), dtype=np.int64)); r = rank_mod(M, p); seen += 1
            best = max(best, r)
            if r == n: break
    else:
        rng = random.Random(p)
        for _ in range(cnt):
            M = sum((rng.randrange(p) * B[j] for j in range(k)), np.zeros((n, n), dtype=np.int64)); r = rank_mod(M, p); seen += 1
            best = max(best, r)
            if r == n: break
    print("p", p, "elements tried", seen, "max rank mod p", best, "(invertible exists)" if best == n else "")
