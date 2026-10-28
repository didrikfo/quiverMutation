"""w+4 only (no w+1..w+3 classification), for width 11 where the full run exceeds 10 min. Usage: theorist_w4only.py WIDTH
Usage: theorist_w4.py WIDTH [WIDTH_B_LIMIT]"""
import sys, collections, time
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
w0 = int(sys.argv[1]); t = time.time()
def tally(d, lens):
    c = collections.Counter()
    for L in lens:
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if lm.matchesAt(L, r, d, ws):
                    o = replay(L, r, d, ws); c[o if isinstance(o, str) else o[1] + '@step%d' % o[0]] += 1
    return c
def kind(c):
    bad = sum(v for k, v in c.items() if k != 'agree')
    return 'none' if not c else ('A' if not bad else ('B' if not c['agree'] else 'C'))
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == w0]
tab = collections.Counter(); tot = collections.Counter()
for d in flo:
    k1 = "skip"; c4 = tally(d, [w0 + 4]); k4 = kind(c4)
    tab[(k1, k4)] += 1; tot.update(c4)
    if k1 != k4: print('CHANGE', d, k1, k4, dict(c4))
print('width', w0, 'rules', len(flo), '(kind w+1..w+3, kind w+4):', dict(tab), 'w+4 matches', dict(tot), 'sec', int(time.time() - t))
