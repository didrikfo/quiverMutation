"""S-1 null: same (i,j) deletion statistic over ALL ordered LNA pairs at n vs same-class pairs. usage: python workshop/rounds/027/maverick_null.py N"""
import sys, os, collections
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh
n = int(sys.argv[1]); lab, _ = classes(n); lab2, _ = classes(n - 1)
ids = {}; M = []; L = []
for m, l in sorted(lab.items()):
    M.append([ids.setdefault(lab2[tuple(pwh.removeVertex(n, list(m), v)[1])], len(ids)) for v in range(1, n + 1)]); L.append(l)
M = np.array(M); K = len(ids); N = len(M)
lid = {l: i for i, l in enumerate(set(L))}; L = np.array([lid[l] for l in L])
for name, (i, j) in (('i=j=mid', ((n + 1) // 2 - 1,) * 2), ('i=1,j=1', (0, 0)), ('i=1,j=n', (0, n - 1))):
    ci = np.bincount(M[:, i], minlength=K); cj = np.bincount(M[:, j], minlength=K)
    allp = ((ci * cj).sum() - (M[:, i] == M[:, j]).sum()) / (N * (N - 1))
    same = tot = 0
    for c in set(L):
        S = M[L == c]; k = len(S)
        a = np.bincount(S[:, i], minlength=K); b = np.bincount(S[:, j], minlength=K)
        same += (a * b).sum() - (S[:, i] == S[:, j]).sum(); tot += k * (k - 1)
    print(n, name, 'all pairs %.3f' % allp, 'same-class pairs %.3f' % (same / tot))
