"""Round 073 (theorist), T5: a finite-dimensional algebra with a cycle meeting D(i), mutated by the repo's procedure.

Parent: vertices 1=i, 2=x, 3=y, 4=z; arrows a:1->2, s:2->3, b:3->4, c:4->3; monomial relations a s, s b c, c b.
TRUE paths: dim e_u L e_v by monomial-algebra enumeration (a path is zero iff it contains a relation as a subword);
finite-dimensional, Cartan det 1, so gl.dim < oo.  The repo's procedure (steps 1-7, no reduce) is run at i = 1; its child
is compared with the child that Oppermann's characterisation (step 7: z_0 = rbar p with (r p)/alpha in I for all alpha,
p a TRUE path) predicts.  Cartan matrices / Coxeter polynomials use true paths (own code, monomial case only).
Run (repo root): .venv/bin/python -W ignore workshop/rounds/073/theorist_cycD.py
"""
import sys, itertools
sys.path.insert(0, ".")
import networkx as nx, sympy as sp
from quivermutation import procedure as pr, arrowPaths as ap

def paths_true(arrows, rels, maxlen=12):
    """All nonzero true paths (as tuples of arrow names) of a monomial algebra; arrows: name->(tail, head)."""
    def zero(p):
        return any(any(p[k:k+len(r)] == r for k in range(len(p)-len(r)+1)) for r in rels)
    verts = sorted({v for t in arrows.values() for v in t})
    out = [((), v, v) for v in verts]
    frontier = [((a,), arrows[a][0], arrows[a][1]) for a in arrows]
    L = 1
    while frontier:
        frontier = [f for f in frontier if not zero(f[0])]
        out += frontier
        if L > maxlen: raise RuntimeError("not finite-dimensional?")
        frontier = [(p + (a,), t, arrows[a][1]) for (p, t, h) in frontier for a in arrows if arrows[a][0] == h]
        L += 1
    return out

def cartan(arrows, rels):
    verts = sorted({v for t in arrows.values() for v in t})
    C = sp.zeros(len(verts))
    for p, t, h in paths_true(arrows, rels):
        C[verts.index(t), verts.index(h)] += 1
    return verts, C

def coxpoly(C):
    Phi = -C.T * C.inv()   # Coxeter matrix, up to transpose convention (same charpoly)
    return sp.factor(Phi.charpoly(sp.Symbol('x')).as_expr())

if __name__ == "__main__":
    A = dict(a=(1, 2), s=(2, 3), b=(3, 4), c=(4, 3))
    R = [("a", "s"), ("s", "b", "c"), ("c", "b")]
    v, C = cartan(A, R); print("parent Cartan", v, C.tolist(), "det", C.det(), "cox", coxpoly(C))
    # repo
    Q = nx.MultiDiGraph(); Q.add_nodes_from([1, 2, 3, 4])
    for n, (t, h) in A.items(): Q.add_edge(t, h, key=0, name=n)
    key = {n: (t, h, 0) for n, (t, h) in A.items()}
    rels = [ap.combination([(tuple(key[x] for x in r), 1)]) for r in R]
    print("gate (vertex-simple, repo):", pr.isMutable(Q, rels, 1))
    Q2, R2 = pr.mutateAtVertex(Q, rels, 1)
    print("repo child arrows:", sorted(Q2.edges(keys=True)))
    print("repo child relations:")
    for r in R2: print("   ", r)
    print("repo child has relation of the form rbar.b.c (a path of length 3 through y,z,y)?  see list above")
    # children after reduce (s = al*.rbar eliminated), arrows named al (2->1), rb (1->3), b, c
    Ch = dict(al=(2, 1), rb=(1, 3), b=(3, 4), c=(4, 3))
    for name, rr in [("repo child (no rbar.b.c)", [("al", "rb", "b", "c"), ("c", "b")]),
                     ("predicted (rbar.b.c = 0 added)", [("rb", "b", "c"), ("c", "b")])]:
        v, C = cartan(Ch, rr)
        print(name, "Cartan", C.tolist(), "det", C.det(), "cox", coxpoly(C), "== parent:", sp.simplify(coxpoly(C) - coxpoly(cartan(A, R)[1])) == 0)
    from quivermutation import search
    print("repo coxeterKey of parent / child (None = guard lets it through):",
          search._coxeterKeyOrNone(pr.toPathAlgebra(Q, rels)), search._coxeterKeyOrNone(pr.toPathAlgebra(Q2, R2)))
    print("search descends from parent? noCycles =", not list(nx.simple_cycles(Q)))
