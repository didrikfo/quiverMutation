"""S-1: distribution of image classes under 'cover-min-mid' and fixed middle vertex, per class, with baseline (all LNAs at n).
usage: python workshop/rounds/027/maverick_dist.py N"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh, quipuForms as qf
n = int(sys.argv[1]); lab, _ = classes(n); lab2, _ = classes(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
def cover(m, v): return sum(1 for s, e in sp(m) if s < v < e)
def mid(m):
    best = min(cover(m, v) for v in range(1, n + 1)); c = [v for v in range(1, n + 1) if cover(m, v) == best]
    return min(c, key=lambda v: (abs(2 * v - n - 1), v))
def img(m, v):
    r = pwh.removeVertex(n, list(m), v); return None if r is None else lab2[tuple(r[1])]
short = {}
def nm(l): return short.setdefault(l, 'T%d' % len(short))
base = collections.Counter(nm(lab2[m]) for m in lab2)
print('n-1 classes sizes', dict(base))
members = collections.defaultdict(list)
for m, l in lab.items(): members[l].append(m)
for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
    d1 = collections.Counter(nm(img(m, mid(m))) if img(m, mid(m)) else None for m in mem)
    d2 = collections.Counter(nm(img(m, (n + 1) // 2)) for m in mem)
    print(str(c)[-30:], len(mem), 'cover-min-mid', dict(d1.most_common(4)), '| v=%d' % ((n + 1) // 2), dict(d2.most_common(4)))
