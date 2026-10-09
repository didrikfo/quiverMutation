"""Reconcile 'key groups' (round 054) with F-047's '25 cospectral groups at n = 10' (round 055).
key = Coxeter polynomial (coxeterTables.lnaKeys).  Counts, for each orbit definition, the number of orbits/classes, polynomial groups
and groups holding >= 2 of them.   usage: python workshop/rounds/055/maverick_recon.py N"""
import sys, collections
sys.path.insert(0, "workshop/rounds/029")
from quivermutation import coxeterTables as ct, freeMoves
n = int(sys.argv[1]); keys = ct.lnaKeys(n)
def uf(lnas, orbits, mirror):
    oid = {m: o for o, mem in orbits.items() for m in mem}; par = {o: o for o in orbits}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    if mirror:
        for m in lnas:
            mr = tuple(freeMoves.mirrorRow(n, m))
            if mr in oid: par[f(oid[m])] = f(oid[mr])
    comps = collections.defaultdict(list)
    for m in lnas: comps[f(oid[m])].append(m)
    return list(comps.values())
def report(name, parts):
    g = collections.defaultdict(list)
    for p in parts: g[keys[p[0]]].append(p)
    print(f"n={n} {name:34s} parts {len(parts):4d}  polynomial groups {len(g):3d}  groups with >=2 parts {sum(len(v)>1 for v in g.values()):3d}")
for name, kw, mir in [("tables+free (table rules + free)", dict(free=True), False),
                      ("+edges", dict(free=True, edges=True), False),
                      ("+edges+doubles", dict(free=True, edges=True, doubles=True), False),
                      ("+edges+doubles+mirror (r054)", dict(free=True, edges=True, doubles=True), True)]:
    lnas, orbits = freeMoves.derivedOrbits(n, rules=None, **kw)
    report(name + "", uf(lnas, orbits, mir))

