"""Round 033 (experimentalist): dim J_i = dim ker g_i on walk rows with J != 0 (gate-admitted).
Usage: experimentalist_dimji.py n class budget_sec [max_exp]   (walk as rounds/027/experimentalist_w.py)
Per (algebra, v) with J != 0: out-degree, {i: dim J_i}, sum, max, number of i with J_i != 0, and the dim e_iAe_v of those i."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3]); maxexp = int(_a[4]) if len(_a) > 4 else 0
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; tab = Counter(); dimtab = Counter(); nalg = 0; stop = False; rows = 0; ex = {}
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
        if time.time() - t0 > budget or (maxexp and nalg >= maxexp): stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            rels = procedure.relationsFrom(alg)
            kd, _ = kerdim(alg, v, rels)
            if kd:
                rows += 1; od = len(ap.arrowsOutOf(alg.quiver, v)); J = perI(alg, v)
                key = (od, tuple(sorted(J.values(), reverse=True))); tab[key] += 1
                for i in J: dimtab[len(ap.allPathsBetween(alg.quiver, i, v)) - len(ap.idealBasis(alg.quiver, rels, i, v))] += 1
                ex.setdefault(key, (v, J, alg.rels))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'class', cls, 'classsize', len(classes[base]), 'expanded', nalg, 'seen', len(seen), 'stopped', stop, '%.0fs' % (time.time() - t0))
print('J != 0 rows', rows)
print('key = (out-degree, sorted dim J_i over i with J_i != 0): rows')
for k, x in sorted(tab.items(), key=str): print(k, x)
print('dim (P_iv/I) at the i with J_i != 0 -> count:', dict(dimtab))
for k, (v, J, r) in ex.items(): print('EX', k, 'v', v, J, r)
