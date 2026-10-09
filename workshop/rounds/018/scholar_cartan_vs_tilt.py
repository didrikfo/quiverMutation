"""T5 round 018: is Cartan congruence (Ladkani Prop 3.6 / Lemma 3.5) the same test as tiltingPlus (AI 2.32(b) = Ladkani 2.3(c))?
Cross-tabulate, at every gate-admitted step of a guarded BFS (as scholar_walk.py) and on the E-078 family,
   tp   = tiltingPlus(parent, v)
   cong = Cartan(reduced child) == R C R^T  (child with the same number of vertices; 'dim' if not)
   guard= child Coxeter key == parent class key
Usage: scholar_cartan_vs_tilt.py n [--class I | --all] [--depth D] [--budget-sec S] | --e078"""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
p = argparse.ArgumentParser(); p.add_argument('n', type=int, nargs='?', default=0)
p.add_argument('--class', type=int, dest='cls', default=-1); p.add_argument('--all', action='store_true')
p.add_argument('--depth', type=int, default=99); p.add_argument('--budget-sec', type=float, default=0, dest='budget')
p.add_argument('--e078', action='store_true'); a = p.parse_args()

tab = Counter(); odd = []
def kerdims(alg, k, rels):
    """dim ker of p |-> (p beta)_beta on paths i~>k mod I, for each i != k"""
    quiver = alg.quiver; outs = ap.arrowsOutOf(quiver, k); res = {}
    for i in quiver.nodes:
        if i == k: continue
        P = ap.allPathsBetween(quiver, i, k)
        if not P: res[i] = 0; continue
        dimV = len(P) - len(ap.idealBasis(quiver, rels, i, k)); rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, v in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): row[(b, kk)] = v
            rows.append(row)
        res[i] = dimV - rank(rows)
    return res
def check(alg, v, base=None, tag=None):
    rels = procedure.relationsFrom(alg)
    if not mutation.mutationIsPossibleAtVertex(alg, v): return None
    t = tiltingPlus(alg.quiver, rels, v)
    raw = mutation.quiverMutationAtVertex(alg, v)
    if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)):
        tab[('illegal', t)] += 1; return None
    ch = reduction.reducePathAlgebra(raw)
    verts = sorted(alg.vertices())
    if len(ch.vertices()) != len(verts): cong = 'dim'
    else: cong = bool((rplus(alg, v, verts).dot(cartan(alg)).dot(rplus(alg, v, verts).T) == cartan(ch)).all())
    guard = None if base is None else (search._coxeterKeyOrNone(ch) == base)
    tab[('tilt' if t else 'NOT', 'cong' if cong is True else 'NOcong(%s)' % cong, 'guard' if guard else ('noguard' if guard is False else '-'))] += 1
    if bool(t) != (cong is True): odd.append((tag, alg.rels, v, t, cong))
    if cong is False:
        D = rplus(alg, v, verts).dot(cartan(alg)).dot(rplus(alg, v, verts).T) - cartan(ch)
        kd = kerdims(alg, v, rels); vi = verts.index(v)
        off = {(i, j): int(D[i, j]) for i in range(len(verts)) for j in range(len(verts)) if D[i, j]}
        rowk = all(i == vi and j != vi for (i, j) in off)
        match = rowk and all(off.get((vi, verts.index(i)), 0) == -d_ for i, d_ in kd.items()) and all(verts.index(i) != vi for i in kd)
        tab[('diff-in-row-k-offdiag', rowk, 'equals -dim ker g_i', match)] += 1
    return ch, guard

if a.e078:
    sys.argv = ['x']
    for n in (5, 6, 7):
        for kind in ('long', 'short', 'zero'):
            for pre in range(0, n - 4):
                post = n - 5 - pre
                A = pathAlgebra.PathAlgebra(); v0 = 1; chain = []
                for _ in range(pre): chain.append(v0); v0 += 1
                aa, b, c, d, e = range(v0, v0 + 5); v0 += 5; tail = list(range(v0, v0 + post))
                A.add_vertices_from(chain + [aa, b, c, d, e] + tail)
                for x, y in zip(chain + [aa], chain[1:] + [aa]): A.add_arrow(x, y)
                for x, y in [(aa, b), (aa, c), (b, d), (c, d), (d, e)]: A.add_arrow(x, y)
                for x, y in zip([e] + tail, tail): A.add_arrow(x, y)
                if kind == 'long': A.add_rel([[aa, b, d, e], [aa, c, d, e]])
                elif kind == 'short': A.add_rel([[aa, b, d], [aa, c, d]])
                else: A.add_rel([[aa, b, d, e]]); A.add_rel([[aa, c, d, e]])
                for v in sorted(A.vertices()):
                    if ap.arrowsOutOf(A.quiver, v): check(A, v, None, (n, kind, pre))
else:
    classes = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
    order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
    sel = range(len(order)) if a.all else [a.cls]
    t0 = time.time()
    for ci in sel:
        base = order[ci]; seen = set(); frontier = []
        def add(alg):
            k = fingerprint.canonicalKey(alg)
            if k is not None:
                if k in seen: return
                seen.add(k)
            frontier.append(alg)
        for s in classes[base]: add(s)
        for dd in range(a.depth):
            cur, frontier[:] = frontier[:], []
            if not cur: break
            for alg in cur:
                if a.budget and time.time() - t0 > a.budget: frontier.clear(); break
                if list(nx.simple_cycles(alg.quiver)): continue
                for v in sorted(alg.vertices()):
                    r = check(alg, v, base, (a.n, ci))
                    if r and r[1]: add(r[0])
        print('class', ci, 'algebras', len(seen), '%.0fs' % (time.time() - t0), flush=True)
for k, v in sorted(tab.items(), key=str): print(k, v)
print('DISAGREEMENTS (tilt xor cong):', len(odd))
for o in odd[:8]: print(o)
