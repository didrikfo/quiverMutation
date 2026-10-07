"""Round 038 referee (experimentalist): test C_B = r C_A r^T + H on walk algebras (n, class, budget s).
H_{v i} = dim J_i (perI), all else 0.  B = reducePathAlgebra(quiverMutationAtVertex) as on the walk.
Usage: experimentalist_review_cartan.py n cls budget_sec"""
import sys, time
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3])
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
from collections import Counter
tab = Counter(); bad = []; vsetbad = 0
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        Q = alg.quiver; V = sorted(Q.nodes)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: vsetbad += 1; continue
            CA = np.array(invariants.cartanMatrix(alg), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
            k = V.index(v); r = np.eye(n, dtype=int)
            for j, w in enumerate(V):
                r[k, j] = -(j == k) + sum(1 for (a, b) in Q.edges() if a == v and b == w)
            # edges() on a MultiDiGraph yields one per parallel arrow
            H = np.zeros((n, n), dtype=int)
            for i, jj in perI(alg, v).items(): H[k, V.index(i)] = jj
            ok = (CB == r @ CA @ r.T + H).all(); ok_noH = (CB == r @ CA @ r.T).all()
            tab[(bool(H.any()), bool(ok), bool(ok_noH))] += 1
            if not ok and len(bad) < 5: bad.append((alg, v, CA, CB, r, H))
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('expanded', nalg, 'vertex-set-mismatch', vsetbad, 'secs %.0f' % (time.time()-t0))
print('(H!=0, identity holds, plain congruence holds):count', dict(tab))
for alg, v, CA, CB, r, H in bad: print('FAIL v', v, '\n', CB - r @ CA @ r.T - H)
