"""S-1: free-vertex deletion (no relation touches v) at length n, split by where v sits: 'head' (left of every relation), 'tail', 'gap'
(between two relations), with the free-run length L the vertex belongs to. For each (source class, type): number of deletions and
whether the image class is one class. Also whether a pair of same-class LNAs with free vertices of the same type give the same image class.
usage: python workshop/rounds/027/maverick_freetype.py N"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh
n = int(sys.argv[1]); lab, _ = classes(n); lab2, _ = classes(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
def img(m, v):
    r = pwh.removeVertex(n, list(m), v); return None if r is None else lab2[tuple(r[1])]
rows = collections.defaultdict(collections.Counter)   # (class, type) -> image counter
for m, l in lab.items():
    if isinstance(l, tuple) and l[1] == '?': continue
    S = sp(m)
    if not S: continue
    first = min(s for s, e in S); last = max(e for s, e in S)
    for v in range(1, n + 1):
        if any(s <= v <= e for s, e in S): continue
        if v < first: t = 'head%d' % (first - 1)
        elif v > last: t = 'tail%d' % (n - last)
        else: t = 'gap'
        rows[(l, t)][img(m, v)] += 1
cls = {}
def nm(c): return cls.setdefault(c, 'C%d' % len(cls))
pure = tot = 0
bytype = collections.defaultdict(lambda: [0, 0])
for (l, t), cnt in sorted(rows.items(), key=lambda t: str(t[0])):
    s = sum(cnt.values()); top = cnt.most_common(1)[0][1]
    ty = t if t == 'gap' else t[:4]
    bytype[ty][0] += top; bytype[ty][1] += s
print('type: purity of image class given (source class, type) -> top/total')
for t, (a, b) in bytype.items(): print(' ', t, a, b, '%.3f' % (a / b))
# head/tail depth profile: purity by exact type string
by2 = collections.defaultdict(lambda: [0, 0])
for (l, t), cnt in rows.items(): by2[t][0] += cnt.most_common(1)[0][1]; by2[t][1] += sum(cnt.values())
print('by exact type', {t: '%d/%d' % tuple(v) for t, v in sorted(by2.items())})
