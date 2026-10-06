"""Round 037 (skeptic): orbits of the E-127 rows (needs skeptic_rows.pkl from skeptic_dump.py) under vertex relabelling and opposite;
import numpy as np
per-row Gamma_i data (d_i, rank of g_i, support of images) against parallel multiplicity and #J_i != 0.
Usage: skeptic_orbits.py skeptic_rows.pkl"""
import sys, pickle, itertools

sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
rows = pickle.load(open(_a[1], 'rb'))
def build(edges, rels):
    Q = nx.MultiDiGraph(); Q.add_nodes_from(range(1, 9))
    for t, h, k in edges: Q.add_edge(t, h, key=k)
    return procedure.toPathAlgebra(Q, [dict(r) for r in rels])
def relabel(alg, pi):
    Q = nx.MultiDiGraph(); Q.add_nodes_from(sorted(pi.values()))
    for t, h, k in alg.quiver.edges(keys=True): Q.add_edge(pi[t], pi[h], key=k)
    m = lambda a: (pi[a[0]], pi[a[1]], a[2])
    rels = [{tuple(m(a) for a in p): c for p, c in r.items()} for r in procedure.relationsFrom(alg)]
    return procedure.toPathAlgebra(Q, rels)
def simple(alg):
    G = nx.DiGraph(); G.add_nodes_from(alg.quiver.nodes)
    for t, h in set((t, h) for t, h, k in alg.quiver.edges(keys=True)): G.add_edge(t, h, m=alg.quiver.number_of_edges(t, h))
    return G
CAP = 10**8
def ck(a):
    k = fingerprint.canonicalKey(a, cap=CAP)
    assert k is not None, 'key None even with cap'
    return k
def isoKeys(alg):
    """All canonicalKeys of alg under quiver-isomorphisms to a fixed labelling = set; use min via iso to itself only for comparison."""
    return None
def same(A, B):
    if int(np.array(A.cartanMatrix()).sum()) != int(np.array(B.cartanMatrix()).sum()) or A.quiver.number_of_edges() != B.quiver.number_of_edges(): return None
    GA, GB = simple(A), simple(B)
    kB = ck(B)
    for m in DiGraphMatcher(GA, GB, edge_match=lambda a, b: a['m'] == b['m']).isomorphisms_iter():
        if ck(relabel(A, m)) == kB: return m
    return None
algs = []
for (nalg, v, edges, rels, J) in rows:
    A = build(edges, rels); algs.append((nalg, v, A, J))
    print('row', nalg, 'v', v, 'J', J, 'nedges', len(edges), 'nrels', len(rels), 'key None?', fingerprint.canonicalKey(A, cap=CAP) is None)
print('--- pairwise: identical up to relabelling (I) / opposite (O)')
for a, b in itertools.combinations(range(len(algs)), 2):
    A, B = algs[a][2], algs[b][2]
    i = same(A, B); o = same(pathAlgebra.dualPathAlgebra(A), B)
    cA, cB = sorted(map(tuple, A.cartanMatrix().tolist())) if hasattr(A.cartanMatrix(), 'tolist') else None, None
    print(algs[a][0], algs[b][0], 'iso' if i else '-', 'opp' if o else '-', 'v-map' , i)
print('--- Cartan e_iAe_v vs d_i (independent check)')
for nalg, v, A, J in algs:
    C = np.array(A.cartanMatrix(), dtype=int)
    print(nalg, 'v', v, {i: (int(C[i-1][v-1]), int(C[v-1][i-1])) for i in J}, 'quiver', sorted(set((t,h,sum(1 for _ in A.quiver.edges(t,h))) for t,h in A.quiver.edges())))
print('--- Cartan multiset invariants')

for nalg, v, A, J in algs:
    C = np.array(A.cartanMatrix(), dtype=int); print(nalg, 'dim', C.sum(), 'det', round(np.linalg.det(C)), 'sorted rowsums', sorted(C.sum(1)), 'sorted colsums', sorted(C.sum(0)))
print('--- Gamma data')
for nalg, v, A, J in algs:
    Q = A.quiver; rels = procedure.relationsFrom(A); outs = ap.arrowsOutOf(Q, v)
    from collections import Counter
    mult = Counter(b[1] for b in outs)
    print('row', nalg, 'v', v, 'outs', outs, 'parallel mult', dict(mult), '#J_i!=0', len(J), 'dimJ', sum(J.values()))
    for i in sorted(Q.nodes):
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        dimV = len(P) - len(ap.idealBasis(Q, rels, i, v)); rws = []
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1])).items(): row[(b, kk)] = x
            rws.append(row)
        tgt = {}
        for t in outs:
            tgt[t] = len(ap.allPathsBetween(Q, i, t[1])) - len(ap.idealBasis(Q, rels, i, t[1])) if i != t[1] else 1
        rk = rank(rws) if rws else 0
        # nonzero images: number of terms per path
        print('   i', i, 'd_i', dimV, 'rank g_i', rk, 'dimJ_i', J.get(i, 0), 'terms per path', sorted(len(r) for r in rws), 'dim e_iAe_t per out', [tgt[t] for t in outs])
