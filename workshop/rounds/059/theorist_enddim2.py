"""Round 059 (theorist): depth-2 version of theorist_enddim.py.  Parent = reduced child of unguarded steps (depth d) from an LNA (any vertex, J = 0 or not), so parents have
parallel arrows and are not derived equivalent to the start in general; the lemma is per step on an arbitrary parent.  Only acyclic parents (the rewrite is not implemented for cycles).
Compares labelled Cartan(child) with dim Hom_K(T_i,T_j) (independent complex code, round 050) and R C R^T - Cartan(child) with -dim Hom(T,T[-1]).
Usage: .venv/bin/python workshop/rounds/059/theorist_enddim2.py N depth cap   (frontier capped at cap per level, seed 59; parallel-arrow parents collected at every level)"""
import sys
sys.path.insert(0, '.')
from collections import Counter
import networkx as nx
_st = open('workshop/rounds/050/skeptic_tilt.py').read()
exec(_st[_st.index('from fractions'):_st.index('def stepTest')])
exec(_st[_st.index('def stepTest'):_st.index('def stepTest') + 400].split('\n\nclass')[0])
from quivermutation import nakayama as nk, pathAlgebra, reduction
N = int(sys.argv[1]); DEPTH = int(sys.argv[2]); cap = int(sys.argv[3]); import random; random.seed(59)
cnt = Counter(); bad = []; parents = 0; seen = set()
def parallel(q):
    return any(q.number_of_edges(a, b) > 1 for a, b in set(q.edges()))
def key(a):
    return (tuple(sorted(a.quiver.edges(keys=True))), tuple(sorted(str(sorted(r.items())) for r in procedure.relationsFrom(a))))
def acyclic(q):
    return not any(True for _ in nx.simple_cycles(nx.DiGraph(q)))
def children(par):
    q = par.quiver; rels = procedure.relationsFrom(par)
    for v in sorted(q.nodes):
        if not ap.arrowsOutOf(q, v): continue
        try:
            nq, nr = procedure.mutateAtVertex(q, rels, v)
            ch = reduction.reducePathAlgebra(procedure.toPathAlgebra(nq, nr))
        except Exception: cnt['err0'] += 1; continue
        yield v, ch, nq, nr
frontier = []
for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
    frontier += [lna, pathAlgebra.dualPathAlgebra(lna)]
allparents = []
for d in range(DEPTH):
    nxt = []
    for par in frontier:
        for v, ch, nq, nr in children(par):
            if not acyclic(ch.quiver): cnt['cyclic child skipped'] += 1; continue
            k = key(ch)
            if k in seen: continue
            seen.add(k); nxt.append(ch)
    random.shuffle(nxt); frontier = nxt[:cap]
    allparents += [a for a in frontier if parallel(a.quiver)]
print('depth', DEPTH, 'frontier', len(frontier), 'parents with parallel arrows', len(allparents))
for par in allparents:
    q = par.quiver; rels = procedure.relationsFrom(par); V = sorted(q.nodes); n = len(V)
    for v in V:
        if not ap.arrowsOutOf(q, v): continue
        try:
            cq, cr = procedure.mutateAtVertex(q, rels, v)
            C = ap.cartanMatrix(cq, cr, exact=True)
            D = procedure.cartanDiscrepancy(q, rels, v, cq, cr)
            _, res = stepTest(par, v, -1)
        except Exception as e:
            cnt['err1'] += 1; continue
        H, Hm, Hp = res[0], res[-1], res[1]
        eq = all(C[i][j] == H[i][j] for i in range(n) for j in range(n)) or all(C[i][j] == H[j][i] for i in range(n) for j in range(n))
        J0 = sum(map(sum, Hm)) == 0 and sum(map(sum, Hp)) == 0
        dm = all(D[i][j] == -Hm[i][j] for i in range(n) for j in range(n)) or all(D[i][j] == -Hm[j][i] for i in range(n) for j in range(n))
        cnt[('J0' if J0 else 'Jne0', 'Hp1=0' if sum(map(sum, Hp)) == 0 else 'Hp1!=0', 'End=child' if eq else 'MISMATCH', 'D=-Hm' if dm else 'D!=-Hm', 'outdeg>1' if len(ap.arrowsOutOf(q, v)) > 1 else 'outdeg1')] += 1
        if not eq: bad.append((v,))
for k, c in sorted(cnt.items(), key=str): print(c, k)
print('mismatches', len(bad))
