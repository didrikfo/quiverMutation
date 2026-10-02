"""Round 023 (skeptic): is a hand-built off-walk rejecting parent reachable by the guarded walk of rounds/022 (same BFS, same
guard, class of the parent's Coxeter key, canonical-key membership)? Usage: skeptic_reach.py n budget_sec
Parents are the list PARENTS below (edges, relations, v) found by skeptic_offwalk.py."""
import sys, time
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
n, budget = int(_argv[1]), float(_argv[2])
PARENTS = {
 6: [("2-out, abde=acde & abdf=acdf", [(1,2),(1,3),(2,4),(2,5),(3,4),(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]],[[1,2,4,6],[1,3,4,6]]]),
     ("shared suffix abdef=acdef", [(1,2),(1,3),(1,4),(2,4),(3,4),(3,6),(4,5),(4,6),(5,6)], [[[1,2,4,5,6],[1,3,4,5,6]]]),
     ("E-078 n=6 (control, 1-out)", [(1,2),(1,3),(2,4),(3,4),(4,5),(5,6)], [[[1,2,4,5],[1,3,4,5]]])],
}
def build(edges, rels):
    A = pathAlgebra.PathAlgebra(); A.add_vertices_from(range(1, n + 1))
    for e in edges: A.add_arrow(*e)
    for r in rels: A.add_rel(r)
    return A
targets = {}
for name, edges, rels in PARENTS[n]:
    A = build(edges, rels); targets[fingerprint.canonicalKey(A)] = (name, search._coxeterKeyOrNone(A))
print({v[0]: str(v[1])[:60] for v in targets.values()})
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
for k0, (nm, ck) in targets.items():
    print(nm, 'class index', order.index(ck) if ck in order else 'NOT A CLASS OF AN LNA/dual', 'size', len(classes.get(ck, [])))
for ci, base in enumerate(order):
    if not any(ck == base for (_, ck) in targets.values()): continue
    t0 = time.time(); seen = set(); frontier = []; found = set(); nexp = 0
    def add(alg):
        k = fingerprint.canonicalKey(alg)
        if k is not None:
            if k in seen: return
            seen.add(k)
            if k in targets: found.add(targets[k][0])
        frontier.append(alg)
    for s in classes[base]: add(s)
    stop = False
    while frontier and not stop:
        cur, frontier[:] = frontier[:], []
        for alg in cur:
            if time.time() - t0 > budget: stop = True; break
            nexp += 1
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == base: add(ch)
    print('class', ci, 'expansions', nexp, 'algebras', len(seen), 'complete' if not stop else 'STOPPED', 'found:', sorted(found), '%.0fs' % (time.time() - t0))
