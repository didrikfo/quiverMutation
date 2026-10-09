"""Round 041 (theorist): random acyclic algebras (not LNA-related). For each gate-admitted, legal step with J != 0 (any support),
record whether the child's Coxeter key equals the PARENT's key (the parent need not be on a walk).  Tests the statement
"J != 0 and gate-admitted  =>  key(child) != key(parent)" outside the LNA class.
Usage: theorist_random.py n budget_sec seed   (run from repo root)"""
import sys, time, random, itertools
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, budget, seed = int(_a[1]), float(_a[2]), int(_a[3]); random.seed(seed)
from collections import Counter
LK = set()
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): LK.add(search._coxeterKeyOrNone(alg))
tab = Counter(); keep = []; t0 = time.time(); tried = 0
def gen():
    nodes = list(range(1, n + 1))
    p = random.choice([0.3, 0.4, 0.5])
    arrows = [(a, b) for a in nodes for b in nodes if a < b and random.random() < p and (b - a <= 3)]
    G = nx.DiGraph(arrows); G.add_nodes_from(nodes)
    allp = {}
    for a in nodes:
        for b in nodes:
            if a < b: allp[(a, b)] = [pth for pth in nx.all_simple_paths(G, a, b)] if a in G and b in G else []
    rels = []
    for (a, b), P in allp.items():
        P2 = [q for q in P]
        if len(P2) >= 2 and random.random() < 0.7:
            q1, q2 = random.sample(P2, 2)
            if len(q1) >= 2 and len(q2) >= 2 and (len(q1) > 2 or len(q2) > 2): rels.append([q1, q2])
        for q in P2:
            if len(q) >= 3 and random.random() < 0.08: rels.append([q])
    return arrows, rels, nodes
while time.time() - t0 < budget:
    arrows, rels, nodes = gen()
    if not arrows: continue
    try:
        A = build(arrows, rels, nodes)
        if list(nx.simple_cycles(A.quiver)): continue
        if any(ap.isIllegalRelation(A.quiver, rr) for rr in procedure.relationsFrom(A)): continue
        kA = search._coxeterKeyOrNone(A)
        tried += 1
        for v in nodes:
            if not mutation.mutationIsPossibleAtVertex(A, v): continue
            J = perI(A, v)
            raw = mutation.quiverMutationAtVertex(A, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != nodes: tab['vertex set changed'] += 1; continue
            same = (search._coxeterKeyOrNone(ch) == kA)
            sh = (len(J), tuple(sorted(J.values())))
            tab[(sh, 'same key' if same else 'key differs', 'LNAkey' if kA in LK else 'noLNAkey')] += 1
            if same and J and kA in LK: keep.append((arrows, rels, v, J))
    except Exception as e:
        tab['err ' + type(e).__name__] += 1
print('n', n, 'seed', seed, 'algebras tried', tried, 'secs %.0f' % (time.time() - t0))
for k, c in sorted(tab.items(), key=str): print('  ', k, c)
print('J != 0, SAME key, key is an LNA key:', len(keep))
for k in keep[:5]: print(k)
