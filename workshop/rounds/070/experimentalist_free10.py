"""S-1 item 1, free-stretch rule (E-114/E-194/E-200) at n = 10 -> 9.
Source label at n = 10: maverick_classes.classes(10) (Coxeter key, cospectral keys split by free+edge+double orbit and mirror; 176 LNAs unresolved, dropped).
Image class at n = 9: exact classifyLength labels (round-065 cached CSV, 20 classes, all 1430 rows).
For each deleted vertex i at distance d >= D from every relation support: per source class, the set of image classes.
usage: .venv/bin/python workshop/rounds/070/experimentalist_free10.py   (about 40 s)"""
import sys, collections, csv
sys.path.insert(0, 'workshop/rounds/027')
from maverick_classes import classes
from quivermutation import nakayama, piecewiseHereditary as pwh
N = 10
lab, bad = classes(N)
hi9 = {r['relations']: r['class'] for r in csv.DictReader(open('workshop/rounds/065/experimentalist_classes_n9.csv'))}
def dist(rl, i):
    ds = [0 if v <= i <= v + r else min(abs(i - v), abs(i - v - r)) for v, r in enumerate(rl, 1) if r]
    return min(ds) if ds else 99
def kind(rl, i):
    sp = [(v, v + r) for v, r in enumerate(rl, 1) if r]
    if not sp: return 'none'
    if i < min(s for s, e in sp): return 'left'
    if i > max(e for s, e in sp): return 'right'
    return 'gap'
badset = set(bad)
members = collections.defaultdict(list)
for m, l in lab.items():
    if m not in badset: members[l].append(m)
imgcache = {}
def img(m, i):
    k = (m, i)
    if k not in imgcache:
        r = pwh.removeVertex(N, list(m), i)
        imgcache[k] = hi9[nakayama.LinearNakayamaAlgebra(r[0], r[1]).relationString()]
    return imgcache[k]
print('n=10 LNAs', len(lab), 'unresolved dropped', len(bad), 'classes', len(members), 'sizes max', max(map(len, members.values())))
res = {}
for D in (1, 2, 3, 4):
    ok = nok = 0; cells = 0; nofit = []; fails = []; gapcells = 0; gapfail = []
    for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
        S = set(); k = 0; Sg = set(); kg = 0
        for m in mem:
            for i in range(1, N + 1):
                if dist(m, i) >= D:
                    t = img(m, i); S.add(t); k += 1
                    if kind(m, i) == 'gap': Sg.add(t); kg += 1
        cells += k; gapcells += kg
        if k == 0: nofit.append((str(c)[-26:], len(mem)))
        elif len(S) == 1: ok += 1
        else: nok += 1; fails.append((str(c)[-26:], len(mem), k, len(S)))
        if kg and len(Sg) > 1: gapfail.append((str(c)[-26:], kg, len(Sg)))
    print('\nD >= %d: classes with a cell %d, one image class %d, not %d; cells %d (gap cells %d); classes with no such cell %d' % (D, ok + nok, ok, nok, cells, gapcells, len(nofit)))
    if fails: print('  failing (class, rows, cells, #images):', fails)
    if gapfail: print('  gap-only cells with >1 image class:', gapfail)
    if D == 3:
        print('  no D>=3 cell (unfitted by the rule): %s' % nofit)
        # same-class pairs, all D>=3 choices agree
        agree = tot = 0
        for c, mem in members.items():
            per = [frozenset(img(m, i) for i in range(1, N + 1) if dist(m, i) >= 3) for m in mem]
            per = [p for p in per if p]
            for a in range(len(per)):
                for b in range(a + 1, len(per)):
                    tot += 1; agree += (per[a] == per[b] and len(per[a]) == 1)
        print('  same-class pairs with D>=3 vertex on both rows: %d, all choices agree on one image: %d' % (tot, agree))
        # image class histogram per source class
        for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
            h = collections.Counter(img(m, i) for m in mem for i in range(1, N + 1) if dist(m, i) >= 3)
            if h: print('   %-26s rows %5d cells %5d -> %s' % (str(c)[-26:], len(mem), sum(h.values()), dict(h)))
