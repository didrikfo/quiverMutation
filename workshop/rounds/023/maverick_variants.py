"""T6 (round 023, maverick): variants of the peeling formula over ALL LNAs n = 8, 9, 10 (depth from workshop/rounds/022/theorist_blocked_depths.txt,
1 for every LNA not listed there, which D1 licenses: 0 mismatches in E-101). Variant v: the chain on a side continues through a relation r iff pred_v(r).
usage: maverick_variants.py   (repository root)"""
import sys
from quivermutation import coxeterTables as ct
data = {}
for line in open('workshop/rounds/022/theorist_blocked_depths.txt'):
    n, d, _, dep = line.split(); data[(int(n), d)] = int(dep)
def make(pa, pb):
    def f(d):
        r = [0] + list(map(int, d)) + [0] * 12; best = None
        for s in range(1, len(d) + 1):
            m = r[s]
            if m < 3: continue
            a = 0
            while s - 1 - a >= 1 and pa(r[s - 1 - a]): a += 1
            b = 0
            while pb(r[s + m - 1 + b]): b += 1
            v = 1 + min(a, b); best = v if best is None else min(best, v)
        return best
    return f
V = {'two/two': make(lambda x: x == 2, lambda x: x == 2), 'any/any': make(lambda x: x >= 2, lambda x: x >= 2),
     'two/any': make(lambda x: x == 2, lambda x: x >= 2), 'any/two': make(lambda x: x >= 2, lambda x: x == 2)}
for n in (8, 9, 10):
    lst = ["".join(map(str, r)) for r in sorted(ct.lnaStatus(n))]
    for name, f in V.items():
        for cls in ('one', 'two+'):
            tot = bad = 0; ex = []
            for d in lst:
                nb = sum(1 for c in d if int(c) >= 3)
                if nb == 0 or (nb == 1) != (cls == 'one'): continue
                tot += 1
                if f(d) != data.get((n, d), 1): bad += 1; ex.append(d)
            print(n, name, cls, 'LNAs', tot, 'mismatch', bad, ex[:6])
