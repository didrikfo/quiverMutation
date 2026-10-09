"""Round 038 (experimentalist): rerun the E-133 walk (n = 8, class idx c0, key-preserving BFS, acyclic algebras, gate-admitted v) and
save, for every (alg, v, i) with path i->v and (J_i != 0 or d_i >= 3): d_i, dim J_i, outdeg(i), outdeg(v), dims e_iAe_v, Cartan row of i,
dim A, distinct-algebra id (canonicalKey, else structural hash), BFS level, expansion index.
Usage: experimentalist_d3table.py n cls budget_sec max_exp out.pkl"""
import sys, time, pickle, hashlib
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget, maxexp, out = int(_a[1]), int(_a[2]), float(_a[3]), int(_a[4]), _a[5]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False; level = 0; rows = []; algs = {}; allrows = 0
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
        Q = alg.quiver; rels = procedure.relationsFrom(alg); V = sorted(Q.nodes)
        cart = {}
        for a in V:
            for b in V:
                if a == b: continue
                P = ap.allPathsBetween(Q, a, b)
                cart[(a, b)] = (len(P) - len(ap.idealBasis(Q, rels, a, b))) if P else 0
        dimA = len(V) + sum(cart.values())
        ck = fingerprint.canonicalKey(alg)
        aid = ('K', str(ck)) if ck is not None else ('H', hashlib.md5(repr((sorted(Q.edges(keys=True)), sorted(cart.items()))).encode()).hexdigest()[:10])
        algs.setdefault(aid, dimA)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            J = perI(alg, v)
            for i in V:
                if i == v or not ap.allPathsBetween(Q, i, v): continue
                d = cart[(i, v)]; j = J.get(i, 0); allrows += 1
                if j or d >= 3:
                    rows.append(dict(exp=nalg, level=level, aid=aid, dimA=dimA, v=v, i=i, d=d, j=j, outi=len(ap.arrowsOutOf(Q, i)), outv=len(ap.arrowsOutOf(Q, v)),
                        outi_targets=sorted(b[1] for b in ap.arrowsOutOf(Q, i)), cartrow=[cart.get((i, b), 0) if b != i else 1 for b in V], V=V, cls=cls))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
    level += 1
pickle.dump(dict(rows=rows, algs=algs, nalg=nalg, allrows=allrows, stopped=stop, secs=time.time()-t0, n=n, cls=cls), open(out, 'wb'))
print('n', n, 'cls', cls, 'expanded', nalg, 'levels', level, 'stopped', stop, 'rows saved', len(rows), 'of', allrows, 'distinct algs seen', len(algs), '%.0fs' % (time.time()-t0))
