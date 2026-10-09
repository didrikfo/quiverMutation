"""Round 023 (scholar): 'reject iff long-sided square', directions and hypotheses.
Definitions (k = v the mutated vertex, alpha the arrows out of v):  J = {c in e_iAe_v : c*alpha = 0 in A for all alpha}.
  tiltingPlus False  <=>  J != 0   (g_i not injective, some i).     gate refuses  <=>  J contains the class of a single path.
  'reject' = gate admits and J != 0  (so J has only non-monomial classes).
Part 1 (--hand): hand-built cases for each hypothesis.   Part 2 (--walk n): over gate-admitted steps of a guarded BFS, tabulate
  (J != 0) against (outdeg v == 1) and (long square in alg.rels).  Usage: scholar_longsquare.py --hand | n [--class I] [--budget-sec S]"""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
p = argparse.ArgumentParser(); p.add_argument('n', type=int, nargs='?', default=0)
p.add_argument('--class', type=int, dest='cls', default=0); p.add_argument('--budget-sec', type=float, default=0, dest='budget')
p.add_argument('--hand', action='store_true'); a = p.parse_args()

def kerdim(alg, v, rels):
    """sum_i dim ker g_i and whether some single path is in J."""
    quiver = alg.quiver; outs = ap.arrowsOutOf(quiver, v); tot = 0; mono = False
    for i in quiver.nodes:
        if i == v: continue
        P = ap.allPathsBetween(quiver, i, v)
        if not P: continue
        dimV = len(P) - len(ap.idealBasis(quiver, rels, i, v))
        rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): row[(b, kk)] = x
            rows.append(row)
            if not row and not ap.isInIdeal(quiver, rels, ap.combination([q])): mono = True
        tot += dimV - (rank(rows) if rows else 0)
    return tot, mono
def longSquare(alg, v):
    outs = ap.arrowsOutOf(alg.quiver, v)
    if len(outs) != 1: return False
    e = outs[0][1]
    for rel in alg.rels:
        if len(rel) >= 2 and all(len(q) >= 3 and q[-1] == e and q[-2] == v for q in rel) and len({q[0] for q in rel}) == 1 and len({q[-3] for q in rel}) == len(rel):
            return True
    return False
def build(arrows, rels, nodes):
    A = pathAlgebra.PathAlgebra(); A.add_vertices_from(nodes)
    for x, y in arrows: A.add_arrow(x, y)
    for r in rels: A.add_rel(r)
    return A
tab = Counter(); odd = []
def row(alg, v, tag=''):
    rels = procedure.relationsFrom(alg); gate = mutation.mutationIsPossibleAtVertex(alg, v)
    kd, mono = kerdim(alg, v, rels); t = bool(tiltingPlus(alg.quiver, rels, v))
    return dict(tag=tag, outdeg=len(ap.arrowsOutOf(alg.quiver, v)), gate=gate, kerdim=kd, tiltingPlus=t, longsq=longSquare(alg, v), mono=mono)
if a.hand:
    # a=1,b=2,c=3,v=4,e=5,f=6,x=7(extra path avoiding v)
    sq = [(1,2),(1,3),(2,4),(3,4)]
    cases = {
     'A long square, v->e only':            (sq+[(4,5)], [[[1,2,4,5],[1,3,4,5]]], 4),
     'B short square (relation at v)':      (sq+[(4,5)], [[[1,2,4],[1,3,4]]], 4),
     'C relation one step past e':          (sq+[(4,5),(5,6)], [[[1,2,4,5,6],[1,3,4,5,6]]], 4),
     'D outdeg 2, commutes into both':      (sq+[(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]],[[1,2,4,6],[1,3,4,6]]], 4),
     'E outdeg 2, commutes into one only':  (sq+[(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]]], 4),
     'F monomial (zero) relation a-b-v-e':  (sq+[(4,5)], [[[1,2,4,5]]], 4),
     'G three paths, two through v (A-B-X)':(sq+[(4,5),(1,7),(7,5)], [[[1,2,4,5],[1,7,5]],[[1,3,4,5],[1,7,5]]], 4),
    }
    for name, (arr, rels, v) in cases.items():
        nodes = sorted({x for e in arr for x in e}); A = build(arr, rels, nodes); print(name, row(A, v), 'rels=', A.rels)
else:
    classes = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
    order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
    t0 = time.time(); base = order[a.cls]; seen = set(); frontier = []
    def add(alg):
        k = fingerprint.canonicalKey(alg)
        if k is not None:
            if k in seen: return
            seen.add(k)
        frontier.append(alg)
    for s in classes[base]: add(s)
    while frontier:
        cur, frontier[:] = frontier[:], []
        for alg in cur:
            if a.budget and time.time() - t0 > a.budget: frontier.clear(); break
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                r = row(alg, v)
                tab[('J!=0' if r['kerdim'] else 'J=0', 'tiltingPlus', r['tiltingPlus'], 'outdeg', min(r['outdeg'], 2), 'longsq', r['longsq'], 'mono', r['mono'])] += 1
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == base: add(ch)
    print('n', a.n, 'class', a.cls, 'algebras', len(seen), '%.0fs' % (time.time() - t0))
    for k, x in sorted(tab.items(), key=str): print(k, x)
