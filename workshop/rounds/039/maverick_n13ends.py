"""n = 13 key class of the lone 3 at (4,5) (hits of maverick_n13class.py): free ends K >= 4, image keys at n = 12; move orbits inside the class.
usage: python workshop/rounds/039/maverick_n13ends.py [orbits]"""
import sys, collections, time, glob
from quivermutation import coxeterTables as ct, freeMoves, lnaMoves, edgeMoves, doubleMutation, piecewiseHereditary as pwh
n = 13
rows = []
for f in sorted(glob.glob('workshop/rounds/039/maverick_n13class_*.txt')):
    rows += [tuple(int(c) for c in l.strip()) for l in open(f) if l.strip()]
print('class size', len(rows))
S = set(rows)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
def ends(m):
    R = sp(m)
    if not R: return None
    return min(s for s, e in R) - 1, n - max(e for s, e in R)
tab = collections.Counter(); bykey = collections.defaultdict(list)
for m in rows:
    e = ends(m)
    for side, K, v in (('head', e[0], 1), ('tail', e[1], n)):
        if K >= 3:
            im = tuple(pwh.removeVertex(n, list(m), v)[1])
            k = ct.lnaCoxeterKey(n - 1, im)
            tab[(K if K < 6 else 6, k)] += 1; bykey[k].append((K, side, ''.join(map(str, m))))
ks = sorted({k for _, k in tab}); print('image keys', len(ks))
for K in (3, 4, 5, 6):
    print('K', K if K < 6 else '>=6', {ks.index(k): c for (kk, k), c in tab.items() if kk == K})
for i, k in enumerate(ks): print(i, k, len(bykey[k]), 'e.g.', bykey[k][:2])
K4 = {k for (K, k) in tab if K >= 4}; print('image keys among K>=4 ends:', sorted(ks.index(k) for k in K4))
if len(sys.argv) > 1:
    idx = {m: i for i, m in enumerate(rows)}; par = list(range(len(rows)))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    def un(a, b): a, b = find(a), find(b); par[a] != par[b] and par.__setitem__(a, b)
    t = time.time(); out = 0
    for m in rows:
        cand = list(lnaMoves.rewritesOf(n, list(m), lnaMoves.ALL_MOVES))
        cand.append(freeMoves.stripLengthTwo(m)); cand += [r for r, _ in edgeMoves.rewritesOf(n, m)]
        cand += [r for r, _ in doubleMutation.rewritesOf(n, m)]
        for r in cand:
            r = tuple(r)
            if r in idx: un(idx[m], idx[r])
            else: out += 1
    print('moves done %.0fs; targets outside class %d' % (time.time() - t, out))
    orb = collections.Counter(find(i) for i in range(len(rows))); print('orbits (forward-and-backward union)', sorted(orb.values(), reverse=True)[:10], len(orb))
    lone = {h: find(idx[tuple([0]*h + [3] + [0]*(n-h-1))[:len(rows[0])]]) for h in (4, 5)} if False else None
    for h in (4, 5):
        r = tuple([0]*h + [3] + [0]*(n - 2 - h - 1))
        print('lone 3 h =', h, 'in class', r in idx, 'orbit size', orb[find(idx[r])] if r in idx else None)
    r4 = tuple([0]*4 + [3] + [0]*6); r5 = tuple([0]*5 + [3] + [0]*5)
    print('same orbit:', find(idx[r4]) == find(idx[r5]))
    mir = collections.Counter((find(i), find(idx[tuple(freeMoves.mirrorRow(n, m))])) for i, m in enumerate(rows)); print('orbit/mirror-orbit pairs', len(mir), sorted(mir.values(), reverse=True)[:6])
    tab2 = collections.Counter()
    for i, m in enumerate(rows):
        e = ends(m)
        for side, K, v in (('head', e[0], 1), ('tail', e[1], n)):
            if K >= 4:
                im = tuple(pwh.removeVertex(n, list(m), v)[1])
                tab2[(orb[find(i)], K, ks.index(ct.lnaCoxeterKey(n - 1, im)))] += 1
    print('K>=4 ends by (orbit size, K, image key idx):', dict(sorted(tab2.items())))
    # n = 12 control: lone 3 (4,4) and (3,5) orbit sizes vs K >= 4 ends of key class
