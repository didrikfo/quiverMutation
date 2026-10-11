"""Round 034 (skeptic): out-degree >= 3 gate-admitted rows with J != 0 on the n = 8 c0 walk (E-126 note), mutated with checkCartan=True.
Usage: skeptic_outdeg3.py n class budget_sec max_exp
Per (algebra, v) with out-degree >= 3: J (perI), then mutateAtVertex(checkCartan=True) -> congruent / FAILS; for failures the support of the
Cartan discrepancy R C R^T - C' is compared with the set {i : J_i != 0} (socle reading control). Also tallies J = 0 out-degree >= 3 rows
(controls: Cartan check on the first 200 of them). Same walk as experimentalist_dimji.py."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget, maxexp = int(_a[1]), int(_a[2]), float(_a[3]), int(_a[4])
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False
tab = Counter(); ctrl = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
def cart(alg, v):
    rels = procedure.relationsFrom(alg)
    try: procedure.mutateAtVertex(alg.quiver, rels, v, checkCartan=True); return 'congruent', None
    except procedure.CartanCongruenceError:
        raw = procedure.mutateAtVertex(alg.quiver, rels, v)
        return 'FAILS', raw
for s in classes[base]: add(s)
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget or (maxexp and nalg >= maxexp): stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            outs = ap.arrowsOutOf(alg.quiver, v)
            if len(outs) >= 3:
                rels = procedure.relationsFrom(alg)
                J = perI(alg, v)
                par = len({b[1] for b in outs}) < len(outs)
                if J or ctrl < 200:
                    if not J: ctrl += 1
                    res, _ = cart(alg, v)
                    rq = alg.quiver; 
                    # discrepancy support
                    supp = None
                    if res == 'FAILS':
                        nq, nr = procedure.mutateAtVertex(alg.quiver, rels, v)
                        D = procedure.cartanDiscrepancy(alg.quiver, rels, v, nq, nr)
                        idx = sorted(alg.quiver.nodes)
                        supp = sorted((idx[a], idx[b]) for a in range(len(D)) for b in range(len(D)) if D[a][b])
                    key = (len(outs), 'parallel' if par else 'simple', 'J!=0' if J else 'J=0', res)
                    tab[key] += 1
                    if J: print('ROW exp', nalg, 'v', v, 'outs', outs, 'J', J, 'dimJ', sum(J.values()), res, 'disc support', supp, 'rels', alg.rels, flush=True)
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'class', cls, 'expanded', nalg, 'seen', len(seen), 'stopped', stop, '%.0fs' % (time.time() - t0))
for k, c in sorted(tab.items(), key=str): print(c, k)
