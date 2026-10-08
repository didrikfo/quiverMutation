"""Round 042 (theorist): which Z-level hypothesis forces B_1 = W_iw + W_wi - W_ii = 1?  Random unitriangular integer Z (m x m, entries 0..e, density p),
pairs (w<i) with Y = Z^-1 satisfying a hypothesis; tally B_1 and also the further coefficients of B(x) = adjS_ii - adjS_wi - adjS_iw.
Usage: theorist_zrandom.py m trials [seed]"""
import sys, numpy as np
from collections import Counter
m, T = int(sys.argv[1]), int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
c = Counter()
for _ in range(T):
    Z = np.triu((rng.random((m, m)) < rng.choice([.2, .4, .6])) * rng.integers(1, 3, (m, m)), 1) + np.eye(m, dtype=int)
    Y = np.rint(np.linalg.inv(Z)).astype(int); W = Y.T @ Z @ Y.T
    for w in range(m):
        for i in range(w + 1, m):
            b0 = 1 - Y[i, w] - Y[w, i]; b1 = W[i, w] + W[w, i] - W[i, i]
            c['all pairs: B0, B1', (int(b0), int(b1))] += 0
            if b0 == 0:
                c['B0=0: B1', int(b1)] += 1
                c['B0=0, Z_wi', (int(Z[w, i]), int(b1))] += 1
            if Y[w, i] == 1: c['Y_wi=1: (B0,B1)', (int(b0), int(b1))] += 1
            if Z[w, i] == 1 and Y[w, i] == 1: c['Z_wi=Y_wi=1: (B0,B1)', (int(b0), int(b1))] += 1
for k in sorted(c, key=str):
    if c[k]: print(k, c[k])
