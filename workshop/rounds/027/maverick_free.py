"""S-1: deletion of a vertex that no relation touches ('free'), or that no relation covers in its interior ('uncovered'), at length n.
Per source class: how many members have such a vertex, the image classes for each choice, whether the image class depends only on the class.
Also the pair statistic: same-class pairs (L, L') each with ALL admissible choices among the free vertices: fraction (L,v,L',v') with equal image class.
usage: python workshop/rounds/027/maverick_free.py N"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
from maverick_classes import classes
from quivermutation import piecewiseHereditary as pwh
n = int(sys.argv[1]); lab, _ = classes(n); lab2, _ = classes(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
def touch(m, v): return any(s <= v <= e for s, e in sp(m))
def interior(m, v): return any(s < v < e for s, e in sp(m))
def img(m, v):
    r = pwh.removeVertex(n, list(m), v); return None if r is None else lab2[tuple(r[1])]
members = collections.defaultdict(list)
for m, l in lab.items():
    if not (isinstance(l, tuple) and l[1] == '?'): members[l].append(m)
for name, pred in (('free (no relation touches v)', lambda m, v: not touch(m, v)),
                   ('uncovered (v interior to no relation)', lambda m, v: not interior(m, v))):
    print('==', name)
    T = E = 0; hasFree = 0; N = 0; pureC = 0
    for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
        sets = [[img(m, v) for v in range(1, n + 1) if pred(m, v)] for m in mem]
        withf = [s for s in sets if s]
        N += len(mem); hasFree += len(withf)
        cnt = collections.Counter(t for s in withf for t in s); tot = sum(cnt.values())
        # ordered pairs of distinct members, all choices
        sq = sum(v * v for v in cnt.values()) - sum(sum(collections.Counter(s).values().__iter__().__class__ and [x * x for x in collections.Counter(s).values()]) for s in withf)
        pr = tot * tot - sum(len(s) ** 2 for s in withf)
        T += pr; E += sq
        # class-function test: is there a single image class common to every free choice of every member?
        common = set.intersection(*[set(s) for s in withf]) if withf else set()
        allsame = all(len(set(s)) == 1 for s in withf)
        pureC += len(mem) if (withf and len(withf) == len(mem) and common) else 0
        if len(mem) >= 6: print('  ', str(c)[-28:], 'size', len(mem), 'with such v', len(withf), 'image classes', dict(cnt.most_common(3)), 'common', bool(common), 'each member single image', allsame)
    print('  LNAs with such a vertex %d of %d; pair rate (all choices) %.3f over %d pairs' % (hasFree, N, E / max(T, 1), T))
