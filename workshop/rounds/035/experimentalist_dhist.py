"""Round 035 (experimentalist): histogram of (d_i, dim J_i), d_i = dim e_iAe_v, at every gate-admitted (v, i) with a path i -> v
on a BFS walk from an LNA class (acyclic algebras only, as 034 theorist_dimji.py walk).  Looks for any d_i >= 3 with J_i != 0.
Usage: experimentalist_dhist.py n class budget_sec   (class = index in classes sorted by (size, key), as 033)"""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3])
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
t0 = time.time(); seen = set(); frontier = []; tab = Counter(); depthmax = Counter(); nalg = 0; stop = False; ex = []; level = 0; nrows = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        Q = alg.quiver; rels = procedure.relationsFrom(alg)
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            J = perI(alg, v)
            for i in Q.nodes:
                if i == v: continue
                P = ap.allPathsBetween(Q, i, v)
                if not P: continue
                d = len(P) - len(ap.idealBasis(Q, rels, i, v)); j = J.get(i, 0)
                tab[(d, j)] += 1; depthmax[level] = max(depthmax[level], d)
                if d >= 3 and j: ex.append((level, v, i, d, j, alg.rels))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
    level += 1
print('n', n, 'class', cls, 'classsize', len(classes[base]), 'expanded', nalg, 'seen', len(seen), 'levels', level, 'stopped', stop, '%.0fs' % (time.time() - t0))
print('(d_i, dim J_i): count at gate-admitted (v,i) with path i->v'); print(dict(sorted(tab.items())))
print('max d_i by BFS level:', dict(sorted(depthmax.items())))
print('J != 0 rows with d >= 3:', len(ex)); [print('EX', e) for e in ex[:5]]
print('J != 0 total:', sum(c for (d, j), c in tab.items() if j), ' by d:', {d: sum(c for (dd, j), c in tab.items() if j and dd == d) for d in sorted({k[0] for k in tab})})
