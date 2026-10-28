"""Class-A rules (all matches replay in E-187/E-205 row model at lengths w+1..w+3): replay at lengths w+4..LMAX (rows <= 14)
and report the set of (r_{ws-1}, r_{ws+w}, boundary) contexts read by the replay.  Usage: theorist_Along.py WMAX LMAX"""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
wmax, lmax = int(sys.argv[1]), int(sys.argv[2])
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] <= wmax]
res = collections.Counter(); nA = 0
for d in flo:
    w = d[0]
    base = collections.Counter()
    for L in range(w + 1, w + 4):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if lm.matchesAt(L, r, d, ws): base[replay(L, r, d, ws) if True else 0] += 1
    if not base or set(base) != {'agree'}: continue
    nA += 1; c = collections.Counter()
    for L in range(w + 4, lmax + 1):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if lm.matchesAt(L, r, d, ws):
                    o = replay(L, r, d, ws); c[o if isinstance(o, str) else 'other'] += 1
    res.update(c)
print('class-A rules of width <=', wmax, ':', nA, ' replay at lengths w+4..%d:' % lmax, dict(res))
