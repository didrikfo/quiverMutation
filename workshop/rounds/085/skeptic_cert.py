"""Round 085 (skeptic), T10: non-Cartan inequivalence certificate for 3060000 vs 3304000 (n = 9) from the quipu theorem
(relation-position combinatorics, quipuForms.quipuForAlmostSeparateLNA) carried across lnaMoves orbits, with power table:
for every n = 9 lnaMoves orbit, (a) all theorem-named members agree (soundness of orbit moves vs theorem), (b) key-equal orbit pairs split into
same-name (certified equivalent = positive controls) and different-name (certified inequivalent = F-010 type pairs).
Run from the repository root:  timeout 10m .venv/bin/python -W ignore workshop/rounds/085/skeptic_cert.py [n]"""
import sys, itertools
sys.path.insert(0, '.')
from quivermutation import lnaMoves, nakayama as nk, search
from quivermutation.quipuForms import quipuForAlmostSeparateLNA as qn
n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
rows = [tuple(r) for r in nk.allRelationLengths(n)]
# lnaMoves.orbitOf is a FORWARD closure (moves are not all symmetric), so take connected components of the move graph
par = {}
def f(x):
    while par.setdefault(x, x) != x: par[x] = par[par[x]]; x = par[x]
    return x
for r in rows:
    s0 = ''.join(map(str, r))
    for m in lnaMoves.orbitOf(n, list(r)): par[f(m)] = f(s0)
comp = {}
for r in rows:
    s0 = ''.join(map(str, r)); comp.setdefault(f(s0), []).append(s0)
orbits = list(comp.values()); seen = {m: i for i, o in enumerate(orbits) for m in o}
print('rows', len(rows), 'orbits (components)', len(orbits), 'covered', len(seen))
info = []
for i, mem in enumerate(orbits):
    names = {}
    for m in mem:
        q = qn(n, [int(c) for c in m])
        if q is not None: names.setdefault(q, []).append(m)
    key = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName(mem[0]))
    info.append((names, key))
bad = [i for i, (nm, k) in enumerate(info) if len(nm) > 1]
unnamed = [i for i, (nm, k) in enumerate(info) if not nm]
print('orbits with two different theorem names (soundness failures):', len(bad), ' unnamed orbits:', len(unnamed))
by = {}
for i, (nm, k) in enumerate(info): by.setdefault(k, []).append(i)
same = diff = unk = 0
for k, g in by.items():
    for a, b in itertools.combinations(g, 2):
        na, nb = info[a][0], info[b][0]
        if not na or not nb: unk += 1; continue
        if set(na) == set(nb): same += 1
        else: diff += 1
print('key-equal orbit pairs: same name (equivalent)', same, ' different name (inequivalent)', diff, ' unnamed', unk)
def orb(s): return seen[s]
for s in ['3060000', '3304000', '2400230', '3004000', '3024000', '6000030']:
    if n == 9:
        i = orb(s); print(s, 'orbit', i, 'size', len(orbits[i]), 'names', {k: v[:3] for k, v in info[i][0].items()})
