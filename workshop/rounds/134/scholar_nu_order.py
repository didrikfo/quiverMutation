"""Round 134 scholar: order of nu on K0 (M = Dmat C^-1 from the Cartan data) for a linear Nakayama algebra.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/134/scholar_nu_order.py RELSTRING [RELSTRING ...]"""
import sys, os
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "087"))
from scholar_rankgrowth import Alg
for rel in sys.argv[1:]:
    A = Alg([int(c) for c in rel]); n = A.n
    C = [[Fraction(int(A.hom(i, j))) for j in range(1, n + 1)] for i in range(1, n + 1)]
    D = [[Fraction(int(A.hom(j, i))) for j in range(1, n + 1)] for i in range(1, n + 1)]
    aug = [C[i][:] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = next(r for r in range(c, n) if aug[r][c]); aug[c], aug[pr] = aug[pr], aug[c]
        aug[c] = [x / aug[c][c] for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c]: f = aug[r][c]; aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    M = [[int(sum(D[i][k] * aug[k][n + j] for k in range(n))) for j in range(n)] for i in range(n)]
    I = [[int(i == j) for j in range(n)] for i in range(n)]; X = I; o = neg = None
    for k in range(1, 2001):
        X = [[sum(X[i][t] * M[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        if neg is None and X == [[-v for v in r] for r in I]: neg = k
        if X == I: o = k; break
    print(rel, "n", n, "order of nu on K0", o, "first M^k = -I", neg)
