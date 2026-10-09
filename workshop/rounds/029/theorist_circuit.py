"""Round 029 (theorist): circuit graph Gamma_i of every out-degree-2 gate-admitted row of the n=8 class-0 walk (round-023/026 walk).
Gamma_i: vertices = nonzero classes of e_iAe_{t b} (one copy per out-arrow b); one edge per nonzero class gamma in e_iAe_v joining [gamma b1],[gamma b2]
(zero product -> edge to ground 'g').  Classes = normal forms up to scalar; we check they are linearly independent and each product is a single class.
Per row: kerdim J, J_b1, J_b2 (per i), and for the skeptic's 22 rows (J_b1,J_b2 != 0 at the same i, J = 0) the shape of Gamma_i.
Usage: theorist_circuit.py n budget_sec [class]"""
import sys, time
from collections import Counter
from fractions import Fraction
sys.path.insert(0, '.'); MYARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
sys.argv = MYARGV
n, budget = int(MYARGV[1]), float(MYARGV[2]); cls = int(MYARGV[3]) if len(MYARGV) > 3 else 0
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; rows_ = []
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
        if time.time() - t0 > budget: frontier.clear(); break
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            if len(ap.arrowsOutOf(alg.quiver, v)) == 2: rows_.append((alg, v))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('walk: algebras', len(seen), 'outdeg-2 rows', len(rows_), '%.0fs' % (time.time() - t0), flush=True)

def norm(d):
    if not d: return None
    k0 = min(d, key=str); c = d[k0]
    return tuple(sorted(((str(k), Fraction(x) / Fraction(c)) for k, x in d.items())))
def jdim(Q, rels, i, v, outs, P):
    dimV = len(P) - len(ap.idealBasis(Q, rels, i, v))
    rw = []
    for q in P:
        r = {}
        for bi, b in enumerate(outs):
            for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1])).items(): r[(bi, kk)] = x
        rw.append(r)
    return dimV - (rank(rw) if rw else 0)
def graph(Q, rels, i, v, outs, P):
    """returns edges list [(u1,u2)], flag ok (classes independent, products single classes)"""
    Bv = ap.idealBasis(Q, rels, i, v); cl = {}; ok = True
    for q in P:
        nf = ap.reduceAgainstPivots(ap.combination([q]), Bv)
        if not nf: continue
        key = norm(nf)
        ends = []
        for bi, b in enumerate(outs):
            pr = ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1]))
            if not pr: ends.append('g')
            else:
                if len(pr) > 1: ok = False
                ends.append((bi, norm(pr)))
        if key in cl and cl[key] != tuple(ends): ok = False
        cl[key] = tuple(ends)
    # independence of class vectors
    if cl:
        vecs = []
        for key in cl: vecs.append({k: x for k, x in key})
        if rank(vecs) != len(vecs): ok = False
    return list(cl.values()), ok
def components(edges):
    par = {}
    def f(x):
        par.setdefault(x, x)
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges:
        if a == 'g': a = ('g', id(b))   # separate ground per edge end so we can count; merge below
        # ground is a single vertex 'g'
    par.clear()
    for a, b in edges:
        pa, pb = f(a if a == 'g' else a), f(b if b == 'g' else b); par[pa] = pb
    comps = {}
    for a, b in edges: comps.setdefault(f(a), []).append((a, b))
    return list(comps.values())
def desc(comp):
    V = {x for e in comp for x in e if x != 'g'}; g = any('g' in e for e in comp)
    ng = sum(e.count('g') for e in comp)
    return (len(comp), len(V), ng)  # edges, nonground vertices, ground endpoints

badex = []; dimc = Counter(); tab = Counter(); shapes22 = Counter(); shapesAll = Counter(); bad = 0; examples = []
for alg, v in rows_:
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v)
    kerJ = 0; per = []
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        dimc[len(P) - len(ap.idealBasis(Q, rels, i, v))] += 1
        d = jdim(Q, rels, i, v, outs, P); d1 = jdim(Q, rels, i, v, [outs[0]], P); d2 = jdim(Q, rels, i, v, [outs[1]], P)
        per.append((i, d, d1, d2, P))
        kerJ += d
    for (i, d, d1, d2, P) in per:
        if d1 > 0 and d2 > 0 and d == 0: cat = '22-type'
        elif d > 0: cat = 'J!=0'
        elif d1 > 0 or d2 > 0: cat = 'one-sided'
        else: continue
        E, ok = graph(Q, rels, i, v, outs, P)
        tab[(cat, ok)] += 1
        if not ok:
            bad += 1
            if len(badex) < 3: badex.append((cat, i, v, alg.rels, alg.quiver.number_of_edges()))
            continue
        comps = components([tuple(e) for e in E])
        # circuit graph restricted: components with a circuit = E >= V (+ground) ; report components that are not trees-with-at-most-one-ground
        sig = tuple(sorted(desc(c) for c in comps if desc(c)[0] > 1 or desc(c)[2] > 0 and desc(c)[0] > 1))
        # non-trivial components: edges>1
        sig = tuple(sorted(desc(c) for c in comps if desc(c)[0] > 1))
        shapesAll[(cat, sig)] += 1
        if cat == '22-type':
            shapes22[sig] += 1
            if len(examples) < 4: examples.append((i, v, alg.rels, [sorted(map(str, e)) for e in E][:8], [ap.arrowsOutOf(Q, v)]))
print('per-(row,i) categories (cat, class-basis ok):', dict(tab), 'not ok:', bad)
print('22-type rows (i level): component shapes (edges, nonground vertices, ground endpoints):')
for k, c in shapes22.most_common(): print('  ', k, c)
print('ALL categories x component shapes:')
for k, c in sorted(shapesAll.items(), key=str): print('  ', k, c)
for e in examples[:2]: print('EXAMPLE i=%d v=%d rels=%s' % (e[0], e[1], e[2]))

print('dim e_iAe_v over all (row,i) with a path i->v (outdeg-2 admitted v):', dict(sorted(dimc.items())))
for b in badex: print('NOT-OK example', b)
