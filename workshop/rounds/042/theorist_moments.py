"""Round 042 (theorist): the moment sequence c_k = chi(e_w, F^k e_w) = (Z^-T F^k)_{ww}, k = -4..8, for J != 0 steps with u = e_w - e_i; split by the orbit relation
(type 1: F e_w = e_i; other: print s with F^s e_w = e_i, |s| <= 4).  Also checks c_{-k} = c_{k-1} (Serre) and the formula
   ratio_k = c_k - c_{k+s} - c_{k-s}  for  B/det S = sum_k (-x)^k ratio_k.
Usage: theorist_moments.py dump.pkl"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']; m = n - 1
c = Counter()
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn or len(r['out']) != 1: continue
    V, v = r['V'], r['v']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ r['CA'] @ rm.T; keep = [a for a in range(n) if a != k]
    Z = Cp[np.ix_(keep, keep)]; w = keep.index(V.index(r['out'][0])); i = keep.index(V.index(list(Jn)[0]))
    Y = np.rint(np.linalg.inv(Z)).astype(int); T = Y.T; F = Z @ T; Fi = np.rint(np.linalg.inv(F)).astype(int)
    def P(k):
        M = np.eye(m, dtype=int)
        for _ in range(abs(k)): M = M @ (F if k > 0 else Fi)
        return M
    e = np.eye(m, dtype=int)
    s = [s for s in range(-4, 5) if (P(s)[:, w] == e[i]).all()]
    seq = {kk: int((T @ P(kk))[w, w]) for kk in range(-5, 9)}
    serre = all(seq[-kk] == seq[kk - 1] for kk in range(0, 6))
    c[('s', tuple(s), 'c_-3..c_6', tuple(seq[kk] for kk in range(-3, 7)), 'Serre', serre)] += 1
for k, v in sorted(c.items(), key=str): print(v, k)
