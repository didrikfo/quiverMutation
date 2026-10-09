"""S-1 first sitting: for same-class LNA pairs at length n, how often do delta_i(L), delta_j(L') land in one class at n-1?
usage: python workshop/rounds/027/maverick_delete.py N [--wild A|B]   (repository root; imports maverick_classes.py)
Classes: Coxeter key, cospectral keys split by orbit+mirror, unresolved LNAs ('?') are put with class --wild (default: own class '?'),
 target length classes = key (+ same split)."""
import sys, os, collections
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh
n = int(sys.argv[1]); target = n - 1
lab, bad = classes(n); lab2, bad2 = classes(target)
ids = {}; 
def cid(l): return ids.setdefault(l, len(ids))
lnas = sorted(lab); members = collections.defaultdict(list)
for m in lnas: members[lab[m]].append(m)
img = {}
for m in lnas:
    row = []
    for v in range(1, n + 1):
        r = pwh.removeVertex(n, list(m), v)
        row.append(-1 if r is None else cid(lab2[tuple(r[1])]))
    img[m] = row
N = n
tot = np.zeros((N, N)); eq = np.zeros((N, N)); valid = np.zeros((N, N))
perclass = {}
for c, mem in members.items():
    M = np.array([img[m] for m in mem]); k = len(mem)
    if k < 2: continue
    K = len(ids) + 1
    cnt = np.zeros((N, K))
    for i in range(N):
        for x in M[:, i]:
            if x >= 0: cnt[i, x] += 1
    ce = np.zeros((N, N)); cv = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            same = (cnt[i] * cnt[j]).sum() - ((M[:, i] == M[:, j]) & (M[:, i] >= 0)).sum()
            vi = (M[:, i] >= 0).sum(); vj = (M[:, j] >= 0).sum()
            both = ((M[:, i] >= 0) & (M[:, j] >= 0)).sum()
            ce[i, j] = same; cv[i, j] = vi * vj - both
    pairs = k * (k - 1)
    perclass[c] = (k, ce, cv)
    tot += pairs; eq += ce; valid += cv
np.set_printoptions(linewidth=200, precision=2, suppress=True)
print('n', n, 'classes', len(members), 'unresolved', len(bad), 'sizes', sorted((len(v) for v in members.values()), reverse=True))
print('ordered same-class pairs', int(tot[0, 0]), '; fraction mapping to same n-1 class, rows i (delete in L) cols j (delete in L\'):')
print((eq / np.maximum(valid, 1)))
print('valid fraction (both deletions admissible)'); print(valid / tot)
# per class, best (i,j) and diagonal
print('per class: size, rate at i=j=1..n (diag), best off-diagonal')
for c, (k, ce, cv) in sorted(perclass.items(), key=lambda t: -t[1][0]):
    pr = k * (k - 1); rate = ce / np.maximum(cv, 1)
    d = np.diag(rate); off = np.unravel_index(np.argmax(rate), rate.shape)
    print(str(c)[:40], k, np.round(d, 2).tolist(), 'best', off[0] + 1, off[1] + 1, round(rate[off], 2))
