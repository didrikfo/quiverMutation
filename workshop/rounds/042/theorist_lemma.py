"""Round 042 (theorist): test the sufficient conditions of the Lemma on J != 0 steps with u = e_w - e_i and C'e_v = e_v.
 Lemma (proved in theorist.md): Z unitriangular after a permutation, F = Z Z^-T, Y = Z^-1, F e_w = e_i.  Then
   Q_1 = 0 and Q_2 = 1 - c_1 + c_2 ... precisely  Q_1 = -c_1 ,  Q_2 = 1 + c_2 - c_1 + sigma_1 * Q_1 ,  c_k = (Y^T F^k)_{ww};  c_1 = Y_iw = 0.
 so Q_2 = 1 iff c_2 = sum_b Y_bw (F e_i)_b = 0; sufficient: Y_bw (F e_i)_b = 0 for every b ('termwise').
Prints, over s = 1 records: how many have termwise zero, how many c_2 = 0, how many have F e_i = -e_m.
Usage: theorist_lemma.py dump.pkl"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']; cnt = Counter()
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn or len(r['out']) != 1 or len(Jn) != 1: continue
    V, v = r['V'], r['v']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ r['CA'] @ rm.T; keep = [a for a in range(n) if a != k]
    if not all(Cp[a, k] == 0 for a in keep): continue
    w = keep.index(V.index(r['out'][0])); i = keep.index(V.index(list(Jn)[0])); e = np.eye(n - 1, dtype=int)
    if not (Cp[k, keep] == e[w] - e[i]).all(): continue
    Z = Cp[np.ix_(keep, keep)]; Y = np.rint(np.linalg.inv(Z)).astype(int); F = Z @ Y.T
    if not (F[:, w] == e[i]).all(): cnt['not s=1'] += 1; continue
    fi = F[:, i]; cnt['s=1'] += 1
    cnt['termwise Y_bw (F e_i)_b = 0 for all b', bool((Y[:, w] * fi == 0).all())] += 1
    cnt['c_2 = sum_b Y_bw (F e_i)_b = 0', int(Y[:, w] @ fi)] += 1
    cnt['F e_i = -e_m', bool((fi <= 0).all() and fi.sum() == -1)] += 1
    cnt['support of F e_i', tuple(sorted(int((fi != 0).sum()) for _ in [0]))] += 1
for kx in sorted(cnt, key=str): print(kx, cnt[kx])
