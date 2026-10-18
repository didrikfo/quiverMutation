"""S-1 item 1, round 073: E-190's pair 90000000 vs 00000030 (equivalent at n = 10, not joined by orbit+mirror moves)
under every vertex deletion to n = 9, labelled by the exact round-065 classifyLength(9) cache; D = distance of the
deleted vertex from the nearest relation support (E-206's definition).  Also the mirror pair and the n = 9 -> 8 context.
usage: .venv/bin/python workshop/rounds/073/experimentalist_pair.py   (a few seconds)"""
import csv
from quivermutation import nakayama, piecewiseHereditary as pwh
hi9 = {r['relations']: r['class'] for r in csv.DictReader(open('workshop/rounds/065/experimentalist_classes_n9.csv'))}
N = 10
def dist(rl, i):
    ds = [0 if v <= i <= v + r else min(abs(i - v), abs(i - v - r)) for v, r in enumerate(rl, 1) if r]
    return min(ds) if ds else 99
def mirror(rs):  # relation dual is not the plain reverse; use the repo's helper if available
    return None
res = {}
for s in ('90000000', '00000030'):
    rl = [int(c) for c in s]
    print('\n%s  (n = 10)' % s)
    for i in range(1, N + 1):
        r = pwh.removeVertex(N, rl, i)
        rs = nakayama.LinearNakayamaAlgebra(r[0], r[1]).relationString()
        print('  delete %2d  D=%2s  -> %s  class %s' % (i, dist(rl, i), rs, hi9[rs][-30:]))
        res[(s, i)] = (dist(rl, i), rs, hi9[rs])
for D in (0, 1, 2, 3):
    imgs = {s: {v[2] for (t, i), v in res.items() if t == s and v[0] >= D} for s in ('90000000', '00000030')}
    print('D >= %d: image classes 90000000 %d, 00000030 %d, union %d' % (D, len(imgs['90000000']), len(imgs['00000030']), len(imgs['90000000'] | imgs['00000030'])))
