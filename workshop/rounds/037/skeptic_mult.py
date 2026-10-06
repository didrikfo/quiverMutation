"""Round 037 (skeptic): presentation-independent mirror test for the E-127 rows. T[(i,j,k)] = dim of e_iAe_j * e_jAe_k inside e_iAe_k
(rank of reduced concatenations). Compare A and B under every quiver isomorphism (vertex relabelling), A vs its dual too.
Usage: skeptic_mult.py skeptic_rows.pkl"""
import sys, pickle, itertools
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/037/skeptic_orbits.py').read().split("algs = []")[0]
sys.argv = ARGV; exec(compile(src, 'orb', 'exec'))
def table(A):
    Q = A.quiver; rels = procedure.relationsFrom(A); V = sorted(Q.nodes); T = {}
    B = {}
    for i in V:
        for j in V:
            if i == j: B[(i, j)] = [()]; continue
            P = ap.allPathsBetween(Q, i, j); piv = ap.idealBasis(Q, rels, i, j)
            B[(i, j)] = P
    for i in V:
        for j in V:
            for k in V:
                if len({i, j, k}) < 3 and not (i == j == k): pass
                if not B[(i, j)] or not B[(j, k)] or i == j or j == k: continue
                ib = ap.idealBasis(Q, rels, i, k); rws = []
                for p in B[(i, j)]:
                    for q in B[(j, k)]:
                        rws.append(ap.reduceAgainstPivots(ap.combination([p + q]), ib))
                T[(i, j, k)] = rank(rws) if rws else 0
    return T
def match(TA, GA, TB, GB):
    for m in DiGraphMatcher(GA, GB, edge_match=lambda a, b: a['m'] == b['m']).isomorphisms_iter():
        if all(TB.get((m[i], m[j], m[k]), 0) == v for (i, j, k), v in TA.items()) and sum(TA.values()) == sum(TB.values()): return m
    return None
algs = [build(e, r) for (_, v, e, r, J) in pickle.load(open(ARGV[1], 'rb'))]
nal = [pickle.load(open(ARGV[1], 'rb'))[a][0] for a in range(len(algs))]
Ts = [table(A) for A in algs]; Gs = [simple(A) for A in algs]
Td = [table(pathAlgebra.dualPathAlgebra(A)) for A in algs]; Gd = [simple(pathAlgebra.dualPathAlgebra(A)) for A in algs]
for a, b in itertools.combinations(range(len(algs)), 2):
    print(nal[a], nal[b], 'sum T', sum(Ts[a].values()), sum(Ts[b].values()), 'iso:', match(Ts[a], Gs[a], Ts[b], Gs[b]), 'opp:', match(Td[a], Gd[a], Ts[b], Gs[b]), flush=True)
