"""Round 021 (scholar): E-093's conditional step, and the A5-shape check of E-095.
At every gate-admitted step (guarded BFS as rounds/018/scholar_cartan_vs_tilt.py, or the E-078 family) with child of the same size:
  rowK[i] : child Cartan [k, i]  vs  dim coker g_i        (g_i : p |-> (p beta)_beta, p in e_i A e_k)   -- derived from steps 4, 6
  colK[i] : child Cartan [i, k]  vs  sum_beta dim e_{t beta} A e_i  - dim e_k A e_i   (= dim ker psi_i)  -- derived from step 7
and for each rejecting parent (tiltingPlus False): does it contain the A5 shape at the mutated vertex v
  (v has exactly one outgoing arrow, to e; two distinct b != c with arrows b->v, c->v; a with arrows a->b, a->c)?
Usage: scholar_step7_entries.py n [--class I] [--budget-sec S] | --e078"""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
p = argparse.ArgumentParser(); p.add_argument('n', type=int, nargs='?', default=0)
p.add_argument('--class', type=int, dest='cls', default=0); p.add_argument('--budget-sec', type=float, default=0, dest='budget')
p.add_argument('--e078', action='store_true'); a = p.parse_args()
tab = Counter(); bad = []; shape = Counter(); rejparents = {}; shapeless = []

def dimAB(quiver, rels, i, j):
    if i == j: return 1
    P = ap.allPathsBetween(quiver, i, j)
    return len(P) - len(ap.idealBasis(quiver, rels, i, j)) if P else 0
def kercoker(alg, k, rels):
    quiver = alg.quiver; outs = ap.arrowsOutOf(quiver, k); res = {}
    for i in quiver.nodes:
        if i == k: continue
        P = ap.allPathsBetween(quiver, i, k)
        tot = sum(dimAB(quiver, rels, i, b[1]) for b in outs)
        dimV = (len(P) - len(ap.idealBasis(quiver, rels, i, k))) if P else 0
        rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, v in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): row[(b, kk)] = v
            rows.append(row)
        r = rank(rows) if rows else 0
        res[i] = (dimV - r, tot - r, tot - dimV)   # ker, coker, sum - dim
    return res
def hasA5(quiver, v):
    outs = ap.arrowsOutOf(quiver, v)
    if len(outs) != 1: return False
    ins = [x for x in set(t for (t, h, _) in ap.arrowsInto(quiver, v))]
    for b in ins:
        for c in ins:
            if b < c:
                srcb = {t for (t, h, _) in ap.arrowsInto(quiver, b)}; srcc = {t for (t, h, _) in ap.arrowsInto(quiver, c)}
                if srcb & srcc: return True
    return False
def check(alg, v, base=None, tag=None):
    rels = procedure.relationsFrom(alg)
    if not mutation.mutationIsPossibleAtVertex(alg, v): return None
    t = bool(tiltingPlus(alg.quiver, rels, v))
    raw = mutation.quiverMutationAtVertex(alg, v)
    if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)): return None
    ch = reduction.reducePathAlgebra(raw)
    verts = sorted(alg.vertices())
    guard = None if base is None else (search._coxeterKeyOrNone(ch) == base)
    if len(ch.vertices()) != len(verts): tab['dim-mismatch'] += 1; return ch, guard
    C = cartan(ch); vi = verts.index(v); kc = kercoker(alg, v, rels)
    okrow = all(C[vi, verts.index(i)] == kc[i][1] for i in kc)
    okcol = all(C[verts.index(i), vi] == kc[i][2] for i in kc)
    tab[('tilt' if t else 'NOT', 'rowk=coker', okrow, 'colk=sum-dim', okcol)] += 1
    if not (okrow and okcol): bad.append((tag, alg.rels, v, t))
    if not t:
        pk = fingerprint.canonicalKey(alg) or str(alg.rels)
        if pk not in rejparents:
            rejparents[pk] = hasA5(alg.quiver, v)
            if not rejparents[pk]: shapeless.append((tag, alg.rels, v))
    return ch, guard
if a.e078:
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
                r = check(alg, v, base, (a.n, a.cls))
                if r and r[1]: add(r[0])
    print('class', a.cls, 'algebras', len(seen), '%.0fs' % (time.time() - t0))
for k, v in sorted(tab.items(), key=str): print(k, v)
print('VIOLATIONS of rowk=coker or colk=sum-dim:', len(bad)); [print(o) for o in bad[:6]]
print('distinct rejecting parents:', len(rejparents), ' A5-shaped at v:', sum(rejparents.values()), ' not:', len(shapeless)); [print(o) for o in shapeless[:6]]
