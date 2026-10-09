"""S-1: 'strip a free end vertex' (delete a vertex left of every relation, or right of every relation) from LNAs whose head (or tail) >= K0.
For each source class: the set of image classes over all such LNAs and all such ends (head and tail pooled), by K0 = 1,2,3. A class is
'transported' if that set has one element. usage: python workshop/rounds/027/maverick_endstrip.py N"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh
n = int(sys.argv[1]); lab, _ = classes(n); lab2, _ = classes(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
def img(m, v):
    r = pwh.removeVertex(n, list(m), v); return lab2[tuple(r[1])]
members = collections.defaultdict(list)
for m, l in lab.items():
    if not (isinstance(l, tuple) and l[1] == '?'): members[l].append(m)
ids = {}
for K0 in (1, 2, 3):
    ok = bad = 0; nlna = 0; detail = []; fmap = {}
    for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
        S = set(); k = 0
        for m in mem:
            R = sp(m)
            if not R: continue
            first = min(s for s, e in R); last = max(e for s, e in R)
            if first - 1 >= K0: S.add(img(m, 1)); k += 1
            if n - last >= K0: S.add(img(m, n)); k += 1
        if k == 0: continue
        nlna += k
        if len(S) == 1: ok += 1; fmap[c] = next(iter(S))
        else: bad += 1; detail.append((str(c)[-22:], len(mem), k, len(S)))
    print('K0 =', K0, ': classes with such ends', ok + bad, 'transported (one image class)', ok, 'not', bad, 'ends used', nlna, 'not:', detail)
    if K0 == 3:
        inv = collections.defaultdict(list)
        for c, t in fmap.items(): inv[t].append(c)
        print('   map classes->image classes: injective?', all(len(v) == 1 for v in inv.values()), 'sources', len(fmap), 'targets', len(inv))
