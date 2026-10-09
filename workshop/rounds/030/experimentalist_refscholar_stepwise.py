"""Round 030 REFEREE of scholar.md: stepwise test d(child)<=max(d,M*d), M=outdeg/indeg at v. Derived from experimentalist_dimdepth.py.
Round 030 (experimentalist): max dim e_iAe_v by BFS mutation depth over the capped walk of an n class (E-115 follow-up).
Per algebra (each canonical key once, at its first-found depth): max over ALL ordered pairs (i,v) of dim e_iAe_v (acyclic quiver), and the same
restricted to v of out-degree 2 that admit mutation (the E-115 'rows').  Depth = BFS level from the class's start algebras.
Usage: experimentalist_dimdepth.py n budget_sec class [maxexp]"""
import sys, time
from collections import Counter, defaultdict
sys.path.insert(0, '.'); MYARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec'))
sys.argv = MYARGV
n, budget, cls = int(MYARGV[1]), float(MYARGV[2]), int(MYARGV[3])
maxexp = int(MYARGV[4]) if len(MYARGV) > 4 else 10**9
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; depth = 0; nexp = 0
recs = []
allmax = defaultdict(Counter); rowmax = defaultdict(Counter); first2 = None; first3 = None; nalg = Counter()
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return False
        seen.add(k)
    frontier.append(alg); return True
def measure(alg, d):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); m = 0; mr = 0
    for v in Q.nodes:
        r2 = len(ap.arrowsOutOf(Q, v)) == 2 and mutation.mutationIsPossibleAtVertex(alg, v)
        for i in Q.nodes:
            if i == v: continue
            P = ap.allPathsBetween(Q, i, v)
            if not P: continue
            dm = len(P) - len(ap.idealBasis(Q, rels, i, v))
            m = max(m, dm)
            if r2: mr = max(mr, dm)
    return m, mr
cache={}
def dmax(alg):
    return measure(alg,0)[0]
for s in classes[base]:
    add(s)
stop = False
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget or nexp >= maxexp: stop = True; break
        if list(nx.simple_cycles(alg.quiver)): continue
        nexp += 1
        m, mr = measure(alg, depth); dpar=max(1,m); allmax[depth][m] += 1; rowmax[depth][mr] += 1; nalg[depth] += 1
        if m >= 2 and first2 is None: first2 = (depth, alg.rels)
        if m >= 3 and first3 is None: first3 = (depth, alg.rels)
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base:
                if list(nx.simple_cycles(ch.quiver)): continue
                dch = max(1, dmax(ch)); Q = alg.quiver
                mo = len(ap.arrowsOutOf(Q, v)); mi = len(ap.arrowsInto(Q, v))
                recs.append((depth, dpar, dch, mo, mi, v, alg.rels, ch.rels))
                add(ch)
    depth += 1
print('n=%d class=%d: seen %d, measured %d, last depth %d, complete=%s, %.0fs' % (n, cls, len(seen), nexp, depth, (not frontier and not stop), time.time() - t0))
print('depth: #algebras | max over all pairs (counts of algebras by max) | max over out-2 rows')
for d in sorted(nalg): print(d, nalg[d], dict(sorted(allmax[d].items())), dict(sorted(rowmax[d].items())))
print('first dim>=2:', first2); print('first dim>=3:', first3)

import pickle
pickle.dump(recs, open('workshop/rounds/030/experimentalist_refscholar_recs_c%d.pkl' % cls, 'wb'))
bad = [r for r in recs if r[2] > max(r[1], max(r[3],r[4])*r[1])]
print('edges', len(recs), 'violations (M=max(out,in)):', len(bad))
for r in bad[:5]: print(r)
print('edges with dpar=1 and dch>=2:', sum(1 for r in recs if r[1]==1 and r[2]>=2), ' of which M=1:', sum(1 for r in recs if r[1]==1 and r[2]>=2 and max(r[3],r[4])<=1))
print('max ratio dch/dpar', max(r[2]/r[1] for r in recs), 'edges with M<=1 and dch>dpar:', sum(1 for r in recs if max(r[3],r[4])<=1 and r[2]>r[1]))
