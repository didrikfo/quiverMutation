"""Round 041 (experimentalist): positive-control hunt. Random small acyclic algebras (m vertices, random arrows, random commutativity and zero relations),
every gate-admitted step with J_i != 0 (perI nonempty): is the child's Coxeter key (computed on the reduced child) equal to the parent's?
Usage: experimentalist_control.py m seed budget_sec [pEdge]"""
import sys, time, random, itertools
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
m, seed, budget = int(_a[1]), int(_a[2]), float(_a[3]); pE = float(_a[4]) if len(_a) > 4 else 0.35
rng = random.Random(seed); t0 = time.time(); N = list(range(1, m + 1))
tab = Counter(); hits = []; tried = 0; built = 0
def paths(arr, s, t, L=None):
    out = []; 
    def go(p):
        if p[-1] == t and len(p) > 1: out.append(tuple(p))
        for (a, b) in arr:
            if a == p[-1] and b <= t: go(p + [b])
    go([s]); return out
while time.time() - t0 < budget:
    tried += 1
    arr = [(a, b) for a in N for b in N if a < b and rng.random() < pE]
    rels = []
    for s in N:
        for t in N:
            if s >= t: continue
            P = [p for p in paths(arr, s, t) if len(p) >= 3]
            if len(P) >= 2 and rng.random() < 0.7:
                k = rng.choice([2, len(P)]) if len(P) > 2 else 2
                sel = rng.sample(P, k); rels.append([list(p) for p in sel])
    for (a, b), (c, d) in itertools.combinations(arr, 2):
        pass
    L2 = [(a, b, c) for (a, b) in arr for (b2, c) in arr if b2 == b]
    for p in L2:
        if rng.random() < 0.15: rels.append([list(p)])
    try:
        A = build(arr, rels, N); procedure.relationsFrom(A)
        if any(ap.isIllegalRelation(A.quiver, rr) for rr in procedure.relationsFrom(A)): continue
        pk = search._coxeterKeyOrNone(A)
    except Exception as e:
        continue
    built += 1
    for v in N:
        try:
            if not mutation.mutationIsPossibleAtVertex(A, v): continue
            J = perI(A, v)
            if not J: continue
            raw = mutation.quiverMutationAtVertex(A, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): tab['illegal child'] += 1; continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != N: tab['vertex dropped'] += 1; continue
            ck = search._coxeterKeyOrNone(ch)
        except Exception as e:
            tab['error'] += 1; continue
        tab[('J!=0 step', pk is not None, ck == pk)] += 1
        if ck == pk and pk is not None: hits.append((v, J, pk, arr, rels))
print('m', m, 'seed', seed, 'tried', tried, 'built', built, '%.0fs' % (time.time() - t0)); [print('  ', k, c) for k, c in sorted(tab.items(), key=str)]
for h in hits[:8]: print('HIT', h)
