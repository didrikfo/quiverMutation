"""Are the n = 10 labels that fail at D >= 2 (680 and 262 rows) one exact class?
Step 1: count free+edge+double orbits (joined with mirror) inside each label; an orbit is a set of known equivalences, so one orbit = one class.
Step 2: if several orbits, try meetingPoints between orbit representatives (a meeting is a proof of equivalence).
usage: timeout 10m .venv/bin/python workshop/rounds/070/experimentalist_pool.py"""
import sys, io, contextlib, collections, runpy, time
sys.argv = ['x']
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path('workshop/rounds/070/experimentalist_free10.py')
members, img, dist, N = (g[k] for k in 'members img dist N'.split())
from quivermutation import freeMoves, nakayama, pathAlgebra, search
lnas, orbits = freeMoves.derivedOrbits(N, rules=None, free=True, edges=True, doubles=True)
oid = {m: o for o, mem in orbits.items() for m in mem}
par = {o: o for o in orbits}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for m in lnas:
    mr = tuple(freeMoves.mirrorRow(N, m))
    if mr in oid: par[f(oid[m])] = f(oid[mr])
fails = []
for c, mem in members.items():
    s = {img(m, i) for m in mem for i in range(1, N + 1) if dist(m, i) >= 2}
    if len(s) > 1: fails.append((c, mem))
for c, mem in fails:
    comp = collections.defaultdict(list)
    for m in mem: comp[f(oid[m])].append(m)
    print(str(c)[-30:], 'rows', len(mem), 'orbit+mirror components', len(comp), sorted(map(len, comp.values()), reverse=True))
    reps = [v[0] for v in comp.values()]
    t0 = time.time()
    for r in reps[1:]:
        a = nakayama.LinearNakayamaAlgebra(N, list(reps[0])); b = nakayama.LinearNakayamaAlgebra(N, list(r))
        for d in (2, 3, 4):
            mp = search.meetingPoints(a.pathAlgebra() if hasattr(a, 'pathAlgebra') else a, b.pathAlgebra() if hasattr(b, 'pathAlgebra') else b, d, alsoDual=True, limit=1)
            if mp: break
        print('  rep', r, 'meets rep0 at depth', d if mp else None, round(time.time() - t0))
        if time.time() - t0 > 240: print('  stopped'); break
# Step 3: image set at D >= 2 per component (a component with 2 images fails inside one exact class)
print('\nD >= 2 image sets per orbit+mirror component')
for c, mem in fails:
    comp = collections.defaultdict(list)
    for m in mem: comp[f(oid[m])].append(m)
    for k, v in sorted(comp.items(), key=lambda t: -len(t[1])):
        h = collections.Counter(img(m, i) for m in v for i in range(1, N + 1) if dist(m, i) >= 2)
        print(str(c)[-30:], 'component rows', len(v), dict(h))
# Step 4: the D = 2 split by side/nearest-relation length, inside the largest component of the 680 label only
print('\nD = 2 split by (side, nearest relation lengths), largest component of the 680-row label')
kind = g['kind']
c, mem = fails[0]
comp = collections.defaultdict(list)
for m in mem: comp[f(oid[m])].append(m)
big = max(comp.values(), key=len)
S = collections.defaultdict(collections.Counter)
for m in big:
    for i in range(1, N + 1):
        if dist(m, i) == 2:
            near = [r for v, r in enumerate(m, 1) if r and min(abs(i - v), abs(i - v - r)) == 2]
            S[(kind(m, i), tuple(sorted(set(near))))][img(m, i)] += 1
for k, h in sorted(S.items()): print('   ', k, dict(h))
