"""Round 034 (theorist): why dim J_i = 1 inside a 2-dim e_iAe_v, and Hom(N,N[-1]) on one example.
Lemma (gate): isMutable says every nonzero PATH p into v has p*b != 0 for some out-arrow b, i.e. no path image lies in
J_i = {x in e_iAe_v : x b = 0 for all out-arrows b}.  Paths span e_iAe_v, so J_i != e_iAe_v: dim J_i <= d_i - 1; J_i != 0 forces d_i >= 2; d_i = 2 forces dim J_i = 1.
Usage: theorist_dimji.py walk n class maxexp   |   theorist_dimji.py hand"""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
def dims(alg, v):
    """{i: (d_i, dim J_i)} over i != v with a path to v; also: any path from an out-neighbour t back to v (cycle)?"""
    Q = alg.quiver; rels = procedure.relationsFrom(alg); res = {}
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        d = len(P) - len(ap.idealBasis(Q, rels, i, v)); res[i] = (d, perI(alg, v).get(i, 0))
    back = any(ap.allPathsBetween(Q, b[1], v) for b in ap.arrowsOutOf(Q, v))
    return res, back
if _a[1] == 'hand':
    def build(arrs, rl):
        A = pathAlgebra.PathAlgebra(); vs = sorted({x for a in arrs for x in a}); A.add_vertices_from(vs)
        for x, y in arrs: A.add_arrow(x, y)
        for r in rl: A.add_rel(r)
        return A
    # H1: E-080 long square a b c d e = 1..5, relation abde = acde, v = d = 4
    H1 = build([(1,2),(1,3),(2,4),(3,4),(4,5)], [[[1,2,4,5],[1,3,4,5]]])
    # H2: i=1 -> a_k (2,3,4) -> v=5 -> t1=6,t2=7, p1 b = p2 b = p3 b for both b (chains of 2-term relations)
    ar = [(1,k) for k in (2,3,4)] + [(k,5) for k in (2,3,4)] + [(5,6),(5,7)]
    rl = [[[1,2,5,t],[1,3,5,t]] for t in (6,7)] + [[[1,3,5,t],[1,4,5,t]] for t in (6,7)]
    H2 = build(ar, rl)
    for name, H, v in (("E-080", H1, 4), ("layered m=3", H2, 5)):
        print(name, "gate", mutation.mutationIsPossibleAtVertex(H, v), "(d_i, dim J_i):", dims(H, v)[0], "path back from out-neighbour to v:", dims(H, v)[1])
else:
    n, cls, maxexp = int(_a[2]), int(_a[3]), int(_a[4])
    classes = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
    order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
    seen = set(); frontier = []; nalg = 0; tab = Counter(); viol = 0; cyc = 0
    def add(alg):
        k = fingerprint.canonicalKey(alg)
        if k is not None:
            if k in seen: return
            seen.add(k)
        frontier.append(alg)
    for s in classes[base]: add(s)
    while frontier and nalg < maxexp:
        cur, frontier[:] = frontier[:], []
        for alg in cur:
            if nalg >= maxexp: break
            nalg += 1
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                r, back = dims(alg, v); cyc += back
                for i, (d, j) in r.items():
                    tab[(d, j)] += 1
                    if d >= 1 and j > d - 1: viol += 1
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == base: add(ch)
    print("n", n, "class", cls, "expanded", nalg, "(d_i, dim J_i): count", dict(sorted(tab.items())), "violations of J<=d-1:", viol, "rows with path back from out-neighbour:", cyc)
