"""Round 042 (theorist): check the block structure C' = [[Z,0],[u^T,1]] (v last), Z = C_A on V\\v, and the closed form
Q = -x[(u+j)^T adj(S)(u+j) - u^T adj(S) u], S = xZ+Z^T; the series Q_1, Q_2 = k^T(W+W^T)e_i-type formula; tally the pieces.
Usage: theorist_blocks.py dump.pkl"""
import sys, pickle
import numpy as np
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']
def prep(r):
    V, v, CA = r['V'], r['v'], r['CA']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    return k, rm @ CA @ rm.T
cnt = Counter(); ex = {}
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    k, Cp = prep(r); V = r['V']
    colv = all(Cp[a, k] == 0 for a in range(n) if a != k)
    keep = [a for a in range(n) if a != k]
    Z = Cp[np.ix_(keep, keep)]; ZA = r['CA'][np.ix_(keep, keep)]
    cnt['colv=e_v', colv] += 1
    cnt['Z=CA restr', bool((Z == ZA).all())] += 1
    cnt['detZ', int(round(np.linalg.det(Z)))] += 1
    if not Jn: continue
    cnt['J support size', len(Jn)] += 1
    i = V.index(list(Jn)[0]); ii = keep.index(i)
    u = Cp[k, keep]; kk = -u
    cnt['CB row v - (u+j) zero', bool((r['CB'][k, keep] == u + np.array([Jn.get(V[a], 0) for a in keep])).all())] += 1
    cnt['k>=0', bool((kk >= 0).all())] += 1
    cnt['k_i', int(kk[ii])] += 1
    Y = np.rint(np.linalg.inv(Z)).astype(int)
    W = Y.T @ Z @ Y.T
    s1 = int(kk @ (Y + Y.T)[:, ii]); q2 = int(u @ (W + W.T)[:, ii] + W[ii, ii])
    cnt['s1 = k^T(Y+Y^T)e_i', s1] += 1; cnt['Q2 formula', q2] += 1
    cnt['W_ii', int(W[ii, ii])] += 1; cnt['k^T W e_i', int(kk @ W[:, ii])] += 1; cnt['k^T W^T e_i', int(kk @ W.T[:, ii])] += 1
for kx in sorted(cnt, key=str): print(kx, cnt[kx])
