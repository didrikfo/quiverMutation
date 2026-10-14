"""Cycles and G2 (round 065, theorist).  Run from the repository root:
    timeout 10m .venv/bin/python workshop/rounds/065/theorist_g2cyclic.py

(1) The one MIN>G instance of theorist_g2exact.py (seed 12, i = 3; parent acyclic with parallel
    arrows, gate-admitted, ker f = 0): the child has the 2-cycle 3 <-> 5.  Its Cartan matrix is
    computed with ALL paths (length <= 8; stable), not vertex-simple ones, and compared with R C R^T.
(2) Cyclic Nakayama parents C_n with all relations of length L: what steps 1-7 do (KeyError / gate).
"""
import sys
sys.path.insert(0, ".")
from fractions import Fraction
import networkx as nx
import numpy as np
from quivermutation import procedure as pr, arrowPaths as ap

def allpaths(Q,L):
    arr=ap.arrowsOf(Q); out={}
    cur=[(a,) for a in arr]
    res=list(cur)
    for l in range(2,L+1):
        nxt=[]
        for p in cur:
            for a in arr:
                if a[0]==p[-1][1]: nxt.append(p+(a,))
        cur=nxt; res+=cur
    return res
def trueCartan(Q,rels,L):
    verts=sorted(Q.nodes)
    paths=allpaths(Q,L)
    byse={}
    for p in paths: byse.setdefault((p[0][0],p[-1][1]),[]).append(p)
    for v in verts: byse.setdefault((v,v),[]).insert(0,())
    # ideal elements: u r w with total length <= L
    gens={}
    for r in rels:
        s,t=ap.source(r),ap.target(r)
        pre=[()]+[p for p in paths if p[-1][1]==s]
        post=[()]+[p for p in paths if p[0][0]==t]
        for u in pre:
            for w in post:
                el={}
                for path,c in r.items():
                    q=u+path+w
                    if len(q)<=L: el[q]=el.get(q,0)+c
                    else: el=None;break
                if el:
                    a=(u[0][0] if u else s); b=(w[-1][1] if w else t)
                    gens.setdefault((a,b),[]).append(el)
    C=[[0]*len(verts) for _ in verts]
    top=0
    for (a,b),ps in byse.items():
        rows=[ {k:Fraction(x) for k,x in g.items()} for g in gens.get((a,b),[])]
        # rank
        piv={}
        rk=0
        for row in rows:
            while row:
                h=min(row)
                if h in piv:
                    f=row[h]
                    for k,x in piv[h].items(): row[k]=row.get(k,0)-f*x
                    row={k:x for k,x in row.items() if x}
                else:
                    row={k:x/row[h] for k,x in row.items()}; piv[h]=row; rk+=1; break
        n=len(ps)-rk
        C[verts.index(a)][verts.index(b)]=n
        if any(len(p)==L for p in ps) and n>0:
            # crude check: nonzero at max length -> may not have stabilised
            pass
    return C


def rcr(Q, i, P, n):
    R = np.eye(n, dtype=int)
    R[i - 1] = 0
    R[i - 1][i - 1] = -1
    for (_t, h, _k) in ap.arrowsOutOf(Q, i):
        R[i - 1][h - 1] += 1
    return R @ P @ R.T


E = [(1, 3, 0), (1, 4, 0), (1, 5, 0), (1, 6, 0), (2, 3, 0), (2, 3, 1), (2, 6, 0), (2, 6, 1), (3, 4, 0), (3, 5, 0), (4, 5, 0), (4, 7, 0), (5, 6, 0), (6, 7, 0)]
rv = [{((2, 3, 1), (3, 5, 0), (5, 6, 0)): 1}, {((1, 3, 0), (3, 4, 0)): 1}, {((4, 5, 0), (5, 6, 0), (6, 7, 0)): 1},
      {((3, 4, 0), (4, 5, 0), (5, 6, 0)): 1, ((3, 5, 0), (5, 6, 0)): 2}, {((3, 4, 0), (4, 5, 0)): 1}, {((5, 6, 0), (6, 7, 0)): 1}]
Q = nx.MultiDiGraph()
for t, h, k in E:
    Q.add_edge(t, h, key=k)
rels = [ap.combination(r.items()) for r in rv]
nq, nr = pr.mutateAtVertex(Q, rels, 3)
P = np.array(ap.cartanMatrix(Q, rels, exact=True))
repoDisc = np.array(pr.cartanDiscrepancy(Q, rels, 3, nq, nr))
true8 = np.array(trueCartan(nq, nr, 8)).T
true10 = np.array(trueCartan(nq, nr, 10)).T
print("(1) gate:", pr.isMutable(Q, rels, 3), "| child has cycle:", not nx.is_directed_acyclic_graph(nq))
print("    repo cartanDiscrepancy (vertex-simple paths) nonzero entries:", int(np.count_nonzero(repoDisc)))
print("    R C R^T - trueCartan(all paths, L=8) zero:", not (rcr(Q, 3, P, 7) - true8).any(),
      "| L=10 same:", not (true8 - true10).any())
minimal = pr._minimalGenerators(nq, [dict(r) for r in nr if r])
print("    relations from i* after repo _minimalGenerators:", sum(1 for r in minimal if ap.source(r) == 3), "(Oppermann count g = 1)")

print("(2) cyclic Nakayama parents")
def cyc(n, lens):
    G = nx.MultiDiGraph()
    for v in range(1, n + 1):
        G.add_edge(v, v % n + 1, key=0)
    def walk(s, L):
        p = []; v = s
        for _ in range(L):
            p.append((v, v % n + 1, 0)); v = v % n + 1
        return tuple(p)
    return G, [ap.combination([(walk(s, L), 1)]) for s, L in lens.items() if L >= 2]
tally = {"KeyError": 0, "gate refuses": 0, "gate admits, rewrite runs": 0}
for n, lens in [(3, {1: 4, 2: 4, 3: 4}), (4, {1: 3, 2: 3, 3: 3, 4: 3}), (4, {1: 5, 2: 4, 3: 3, 4: 4}),
                (5, {s: 4 for s in range(1, 6)}), (3, {1: 3, 2: 3, 3: 3})]:
    G, rl = cyc(n, lens)
    for i in range(1, n + 1):
        ok = pr.isMutable(G, rl, i)
        try:
            pr.mutateAtVertex(G, rl, i)
            tally["gate admits, rewrite runs" if ok else "gate refuses"] += 1
        except KeyError:
            tally["KeyError"] += 1
print("   ", tally)
