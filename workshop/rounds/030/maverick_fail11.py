"""Dissect the one K0=3 failure at n=11 of maverick_endstrip2: which source LNAs go to which image class. usage: python workshop/rounds/030/maverick_fail11.py"""
import sys, collections
sys.path.insert(0, "workshop/rounds/030"); sys.path.insert(0, "workshop/rounds/027")
import maverick_endstrip2 as e2, maverick_classes as mc
from quivermutation import piecewiseHereditary as pwh, freeMoves
n = 11; lab, _ = e2.corrected(n); lab2, _ = e2.corrected(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
cls = collections.defaultdict(list)
for m, l in lab.items(): cls[l].append(m)
tgt = [c for c in cls if str(c).endswith('-3, -3, -2, -1, 0, 1, 1)') and len(cls[c]) == 1305][0]
print('class', tgt)
byimg = collections.defaultdict(list)
for m in cls[tgt]:
    R = sp(m); first = min(s for s, e in R); last = max(e for s, e in R)
    for side, cond, v in (('head', first - 1 >= 3, 1), ('tail', n - last >= 3, n)):
        if cond: byimg[lab2[tuple(pwh.removeVertex(n, list(m), v)[1])]].append((side, ''.join(map(str, m))))
for t, L in byimg.items():
    print('image', t, len(L), 'e.g.', L[:4], collections.Counter(s for s, _ in L))
# orbit structure of source class
lnas, orbits = freeMoves.derivedOrbits(n, rules=None, free=True, edges=True, doubles=True)
oid = {m: o for o, mem in orbits.items() for m in mem}
print('source orbits (free+edges+doubles):', collections.Counter(oid[m] for m in cls[tgt]).most_common(6))
tab = collections.Counter()
for m in cls[tgt]:
    R = sp(m); first = min(s for s, e in R); last = max(e for s, e in R)
    for side, cond, v in (('head', first - 1 >= 3, 1), ('tail', n - last >= 3, n)):
        if cond: tab[(oid[m], str(lab2[tuple(pwh.removeVertex(n, list(m), v)[1])])[-30:])] += 1
print(dict(tab))
mir = collections.Counter((oid[m], oid[tuple(freeMoves.mirrorRow(n, m))]) for m in cls[tgt]); print('orbit vs mirror orbit', dict(mir))
print('--- by free-run length K (this class, orbit 15107)')
t2 = collections.Counter()
for m in cls[tgt]:
    if oid[m] != 15107: continue
    R = sp(m); first = min(s for s, e in R); last = max(e for s, e in R)
    for side, K, v in (('head', first - 1, 1), ('tail', n - last, n)):
        if K >= 3: t2[(K, str(lab2[tuple(pwh.removeVertex(n, list(m), v)[1])])[-17:-9])] += 1
for k in sorted(t2): print(k, t2[k])
