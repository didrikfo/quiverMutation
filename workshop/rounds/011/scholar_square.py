"""H-015 / T5: one-map identity on commutative-square-into-a-vertex instances, n = 5..7.

Instance family (vertices relabelled 1..n): optional chain tail t0 > ... > a, square
a>b, a>c, b>d, c>d, then d>e and optional tail e>... ; relations:
  'long'  : abde = acde (relation of length 4, only paths of length 3 survive) -- as E-032 step 7
  'short' : abd  = acd  (true commutative square; then abde = acde follows)
  'zero'  : abde = 0 and acde = 0 (monomial control)
At every vertex k: gate (mutationIsPossibleAtVertex), tiltingPlus (= AI 2.32(b) = Ladkani 2.3(c)),
and, when the gate admits and the child is acyclic, Cartan(child) == R Cartan R^T (independent of
the map: it uses the actual rewrite).
"""
import sys; sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
import networkx as nx

def build(pre, post, kind):
    # vertices: pre chain, then a,b,c,d,e, then post chain
    A = pathAlgebra.PathAlgebra()
    v = 1; chain = []
    for _ in range(pre): chain.append(v); v += 1
    a, b, c, d, e = range(v, v + 5); v += 5
    tail = list(range(v, v + post))
    for x, y in zip(chain + [a], chain[1:] + [a]) if chain else []: A.add_arrow(x, y)
    A.add_vertices_from(chain + [a, b, c, d, e] + tail)
    for x, y in [(a, b), (a, c), (b, d), (c, d), (d, e)]: A.add_arrow(x, y)
    seq = [e] + tail
    for x, y in zip(seq, seq[1:]): A.add_arrow(x, y)
    if kind == 'long': A.add_rel([[a, b, d, e], [a, c, d, e]])
    elif kind == 'short': A.add_rel([[a, b, d], [a, c, d]])
    elif kind == 'zero': A.add_rel([[a, b, d, e]]); A.add_rel([[a, c, d, e]])
    return A, d

print("kind pre post n | k=d: gate tilt | all vertices: (gate,tilt) tally | gate-admitted rejections | admitted not cong")
for n in (5, 6, 7):
    for kind in ('long', 'short', 'zero'):
        for pre in range(0, n - 4):
            post = n - 5 - pre
            alg, d = build(pre, post, kind)
            rels = procedure.relationsFrom(alg)
            verts = sorted(alg.vertices()); tally = Counter(); bad = []
            for k in verts:
                if not ap.arrowsOutOf(alg.quiver, k): continue
                g = mutation.mutationIsPossibleAtVertex(alg, k)
                t = tiltingPlus(alg.quiver, rels, k)
                tally[(g, t)] += 1
                if g:
                    ch = mutation.quiverMutationAtVertex(alg, k)
                    if any(ap.isIllegalRelation(ch.quiver, r) for r in procedure.relationsFrom(ch)): bad.append((k, 'illegal')); continue
                    ch = reduction.reducePathAlgebra(ch)
                    R = rplus(alg, k, verts)
                    cong = bool((R.dot(cartan(alg)).dot(R.T) == cartan(ch)).all())
                    if not (t and cong): bad.append((k, 't=%s cong=%s' % (t, cong)))
            gd = mutation.mutationIsPossibleAtVertex(alg, d); td = tiltingPlus(alg.quiver, rels, d)
            print(kind, pre, post, n, '| d=%d gate=%s tilt=%s' % (d, gd, td), '|', dict(tally), '|', bad)
