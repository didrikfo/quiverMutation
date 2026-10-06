"""Round 037 (skeptic): re-run the E-127 n=8 c0 walk (same as rounds/034/skeptic_outdeg3.py) and pickle the algebras of every
out-degree >= 3 gate-admitted row with J != 0. Usage: skeptic_dump.py budget_sec max_exp out.pkl"""
import sys, time, pickle
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls = 8, 0; budget, maxexp, out = float(_a[1]), int(_a[2]), _a[3]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False; rows = []
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
            outs = ap.arrowsOutOf(alg.quiver, v)
            if len(outs) >= 3:
                J = perI(alg, v)
                if J:
                    rows.append((nalg, v, sorted(alg.quiver.edges(keys=True)), [dict(r) for r in procedure.relationsFrom(alg)], J))
                    print('row', nalg, v, J, flush=True)
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
pickle.dump(rows, open(out, 'wb')); print('done', nalg, len(rows), '%.0fs' % (time.time()-t0))
