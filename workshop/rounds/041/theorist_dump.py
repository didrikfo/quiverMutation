"""Round 041 (theorist): dump (C_A, V, v, J-vector, C_B, passed) for every gate-admitted step on the guarded class walk (as rounds/039/theorist_guard.py).
Usage: theorist_dump.py n cls budget_sec out.pkl   (run from repo root)"""
import sys, time, pickle
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget, out = int(_a[1]), int(_a[2]), float(_a[3]), _a[4]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False; recs = []
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
        Q = alg.quiver; V = sorted(Q.nodes)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: continue
            CA = np.array(invariants.cartanMatrix(alg), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
            J = perI(alg, v)
            passed = search._coxeterKeyOrNone(ch) == base
            recs.append(dict(CA=CA, CB=CB, V=V, v=v, J={i: jj for i, jj in J.items()}, passed=passed,
                             out=[w for (a, w) in Q.edges() if a == v]))
            if passed: add(ch)
pickle.dump(dict(recs=recs, n=n, cls=cls, nalg=nalg), open(out, 'wb'))
print('n', n, 'cls', cls, 'expanded', nalg, 'records', len(recs), 'secs %.0f' % (time.time() - t0))
