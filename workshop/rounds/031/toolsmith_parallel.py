"""Round 031 (toolsmith): parallel-arrow controls for rule W (E-107, E-110, E-111) and an independent check of the
dim e_iAe_v count on quivers with doubled arrows (E-117 caveat).
Usage: toolsmith_parallel.py --controls | --dimcheck [N seed]
Controls are built with arrow-level relations (procedure.toPathAlgebra), so a relation may run along a doubled arrow.
Independent count: enumerate arrow paths by DFS, build the ideal part u*r*w spanning set, exact rank over Q (sympy). No arrowPaths call."""
import sys, random
from fractions import Fraction
sys.path.insert(0, '.'); MYARGS = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = MYARGS
import networkx as nx, sympy

def W(alg, v):   # verbatim from rounds/027/experimentalist_w.py
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v)
    for R in rels:
        if len(R) < 2: continue
        for b1 in outs:
            if not all(len(q) >= 2 and q[-1] == b1 and q[-2][1] == v for q in R): continue
            b2 = [b for b in outs if b != b1][0]
            x = ap.combination({q[:-1]: c for q, c in R.items()})
            if ap.isInIdeal(Q, rels, x): continue
            if ap.isInIdeal(Q, rels, ap.combination({q[:-1] + (b2,): c for q, c in R.items()})): return True
    return False

# ---------- independent count ----------
def paths(Q, i, v):
    out = []
    def go(u, cur):
        if u == v and cur: out.append(tuple(cur))
        for (t, h, k) in sorted(Q.out_edges(u, keys=True)):
            go(h, cur + [(t, h, k)])
    go(i, []); return out
def ideal_rows(Q, rels, i, v, index):
    rows = []
    for R in rels:
        p0 = next(iter(R)); s, t = p0[0][0], p0[-1][1]
        for u in (paths(Q, i, s) if i != s else [()]):
            for w in (paths(Q, t, v) if t != v else [()]):
                row = [0] * len(index)
                for q, c in R.items(): row[index[u + q + w]] += sympy.Rational(c)
                rows.append(row)
    return rows
def dimA(Q, rels, i, v):
    P = paths(Q, i, v)
    if not P: return 0
    idx = {p: n for n, p in enumerate(P)}; rows = ideal_rows(Q, rels, i, v, idx)
    return len(P) - (sympy.Matrix(rows).rank() if rows else 0)
def dimJ(Q, rels, i, v):
    """dim ker(g_i : e_iAe_v -> sum_b e_iAe_{t(b)}), b over arrows out of v."""
    P = paths(Q, i, v); outs = sorted(Q.out_edges(v, keys=True))
    if not P or not outs: return len(P) and dimA(Q, rels, i, v)
    blocks = []; off = 0
    for b in outs:
        Pb = paths(Q, i, b[1]); blocks.append((off, {p: n for n, p in enumerate(Pb)}, Pb)); off += len(Pb)
    tot = off; Irows = []
    for (o, idx, Pb), b in zip(blocks, outs):
        for r in ideal_rows(Q, rels, i, b[1], idx): Irows.append([0] * o + r + [0] * (tot - o - len(Pb)))
    G = []
    for p in P:
        row = [0] * tot
        for (o, idx, Pb), b in zip(blocks, outs): row[o + idx[p + (b,)]] += 1
        G.append(row)
    rI = sympy.Matrix(Irows).rank() if Irows else 0
    rGI = sympy.Matrix(G + Irows).rank()
    return dimA(Q, rels, i, v) - (rGI - rI)

def mk(arrows, rels):
    q = nx.MultiDiGraph()
    for t, h in arrows: q.add_edge(t, h)
    return q, rels
def arrow_ids(arrows):
    """arrows list with repeats -> list of (t,h,k)"""
    cnt = {}; ids = []
    for t, h in arrows: k = cnt.get((t, h), 0); cnt[(t, h)] = k + 1; ids.append((t, h, k))
    return ids

def controls():
    # each: name, arrows, rel builder (given ids by (t,h,k)), v
    def A(t, h, k=0): return (t, h, k)
    C = {}
    # 1. parallel W (positive): doubled 2=>4; commute into b1=4->5; both doubled paths killed by b2=4->6
    C['P1 parallel W: g0 b1=g1 b1, g0 b2=g1 b2=0'] = ([(2,4),(2,4),(4,5),(4,6)],
        [{(A(2,4,0),A(4,5)):1,(A(2,4,1),A(4,5)):-1},{(A(2,4,0),A(4,6)):1},{(A(2,4,1),A(4,6)):1}], 4)
    # 2. same with a prefix arrow 1->2 (dim e1Ae4 = 2 too)
    C['P2 parallel W with prefix 1->2'] = ([(1,2),(2,4),(2,4),(4,5),(4,6)],
        [{(A(1,2),A(2,4,0),A(4,5)):1,(A(1,2),A(2,4,1),A(4,5)):-1},{(A(1,2),A(2,4,0),A(4,6)):1},{(A(1,2),A(2,4,1),A(4,6)):1}], 4)
    # 3. cancels (nn) with doubled arrow: commute into both out-arrows
    C['P3 parallel cancels: g0 b=g1 b for b1 and b2'] = ([(2,4),(2,4),(4,5),(4,6)],
        [{(A(2,4,0),A(4,5)):1,(A(2,4,1),A(4,5)):-1},{(A(2,4,0),A(4,6)):1,(A(2,4,1),A(4,6)):-1}], 4)
    # 4. parallel OUT-arrows at v: square into v, then doubled 4=>5: commute through b0, both killed by b1
    C['P4 parallel out-arrows b0,b1 (4=>5)'] = ([(1,2),(1,3),(2,4),(3,4),(4,5),(4,5)],
        [{(A(1,2),A(2,4),A(4,5,0)):1,(A(1,3),A(3,4),A(4,5,0)):-1},{(A(1,2),A(2,4),A(4,5,1)):1},{(A(1,3),A(3,4),A(4,5,1)):1}], 4)
    # 5. negative: doubled arrow, no relation (everything free): gate? J should be 0
    C['N1 doubled, free'] = ([(2,4),(2,4),(4,5),(4,6)], [], 4)
    # 6. negative: doubled, g0 b1 = g1 b1 only (b2 free): J=0
    C['N2 doubled, commute into b1 only'] = ([(2,4),(2,4),(4,5),(4,6)],
        [{(A(2,4,0),A(4,5)):1,(A(2,4,1),A(4,5)):-1}], 4)
    C['N3 doubled, g0 b1 = 0 (single path)'] = ([(2,4),(2,4),(4,5),(4,6)], [{(A(2,4,0),A(4,5)):1}], 4)
    # 7. out-degree 1 through a doubled arrow (long square): g0 e = g1 e
    C['P5 parallel G: outdeg 1, g0 e = g1 e'] = ([(2,4),(2,4),(4,5)], [{(A(2,4,0),A(4,5)):1,(A(2,4,1),A(4,5)):-1}], 4)
    # 8. tripled arrow, H chain (3-term kernel): g0b1=g1b1, g1b2=g2b2, g0b2=0, g2b1=0
    C['P6 parallel H: tripled arrow chain'] = ([(2,4),(2,4),(2,4),(4,5),(4,6)],
        [{(A(2,4,0),A(4,5)):1,(A(2,4,1),A(4,5)):-1},{(A(2,4,1),A(4,6)):1,(A(2,4,2),A(4,6)):-1},{(A(2,4,0),A(4,6)):1},{(A(2,4,2),A(4,5)):1}], 4)
    # plain (non-parallel) twins, E-110 cases D and W-type, for comparison
    C['D0 plain cancels (E-110 D)'] = ([(1,2),(1,3),(2,4),(3,4),(4,5),(4,6)],
        [{(A(1,2),A(2,4),A(4,5)):1,(A(1,3),A(3,4),A(4,5)):-1},{(A(1,2),A(2,4),A(4,6)):1,(A(1,3),A(3,4),A(4,6)):-1}], 4)
    C['W0 plain W-type'] = ([(1,2),(1,3),(2,4),(3,4),(4,5),(4,6)],
        [{(A(1,2),A(2,4),A(4,5)):1,(A(1,3),A(3,4),A(4,5)):-1},{(A(1,2),A(2,4),A(4,6)):1},{(A(1,3),A(3,4),A(4,6)):1}], 4)
    return C

def run_controls():
    print('%-48s %4s %5s %5s %4s %5s  %s' % ('case','gate','dimA','J','W','dimJ*','Cartan'))
    for name, (arrows, rels, v) in controls().items():
        q = nx.MultiDiGraph()
        for t, h in arrows: q.add_edge(t, h)
        alg = procedure.toPathAlgebra(q, rels); R = procedure.relationsFrom(alg)
        gate = mutation.mutationIsPossibleAtVertex(alg, v)
        kd, mono = kerdim(alg, v, R)
        dims = {i: (len(ap.allPathsBetween(alg.quiver, i, v)) - len(ap.idealBasis(alg.quiver, R, i, v)), dimA(alg.quiver, R, i, v)) for i in alg.quiver.nodes if i != v and ap.allPathsBetween(alg.quiver, i, v)}
        jstar = sum(dimJ(alg.quiver, R, i, v) for i in dims)
        try: procedure.mutateAtVertex(alg.quiver, R, v, checkCartan=True); cart = 'congruent'
        except procedure.CartanCongruenceError: cart = 'FAILS'
        except Exception as e: cart = 'error ' + type(e).__name__
        w = W(alg, v) if len(ap.arrowsOutOf(alg.quiver, v)) == 2 else None
        print('%-48s %4s %5s %5d %4s %5d  %s   (code dim, indep dim) by i: %s mono=%s' % (name, gate, max(d[0] for d in dims.values()), kd, w, jstar, cart, dims, mono))

def dimcheck(N, seed):
    nchild = [0]; cpairs = [0]; cbad = [0]
    rng = random.Random(seed); bad = 0; npairs = 0; npar = 0; nJ = 0; Jpos = 0; Jbad = 0
    for _ in range(N):
        nv = rng.randint(4, 6); arrows = []
        for t in range(1, nv + 1):
            for h in range(t + 1, nv + 1):
                if rng.random() < 0.45: arrows += [(t, h)] * rng.choice([1, 1, 2])
        if not arrows: continue
        q = nx.MultiDiGraph(); q.add_nodes_from(range(1, nv + 1))
        for t, h in arrows: q.add_edge(t, h)
        if not any(q.number_of_edges(t, h) > 1 for t, h in set(arrows)): continue
        npar += 1; rels = []
        for _r in range(rng.randint(1, 4)):
            s = rng.randint(1, nv); t = rng.randint(s + 1, nv) if s < nv else s
            P = [p for p in paths(q, s, t) if len(p) >= 2] if s < t else []
            if not P: continue
            if rng.random() < 0.4 or len(P) == 1: rels.append({rng.choice(P): 1})
            else: p1, p2 = rng.sample(P, 2); rels.append({p1: 1, p2: -1})
        rels = [r for r in rels]
        alg = procedure.toPathAlgebra(q, rels); R = procedure.relationsFrom(alg)
        for i in q.nodes:
            for v in q.nodes:
                if i == v or not ap.allPathsBetween(q, i, v): continue
                a = len(ap.allPathsBetween(q, i, v)) - len(ap.idealBasis(q, R, i, v)); b = dimA(q, R, i, v); npairs += 1
                if a != b: bad += 1; print('DIM MISMATCH', sorted(q.edges(keys=True)), R, i, v, a, b)
        for v in q.nodes:
            if not q.out_edges(v): continue
            kd, _ = kerdim(alg, v, R); kj = sum(dimJ(q, R, i, v) for i in q.nodes if i != v and ap.allPathsBetween(q, i, v))
            nJ += 1; Jpos += kd > 0
            if kd == 0 and mutation.mutationIsPossibleAtVertex(alg, v):
                try: cq, cr = procedure.mutateAtVertex(q, R, v)
                except Exception: continue
                if not nx.is_directed_acyclic_graph(nx.DiGraph(cq)) or any(len(cq.out_edges(x)) and False for x in cq): continue
                nchild[0] += 1
                for i in cq.nodes:
                    for w in cq.nodes:
                        if i == w or not ap.allPathsBetween(cq, i, w): continue
                        a = len(ap.allPathsBetween(cq, i, w)) - len(ap.idealBasis(cq, cr, i, w)); b = dimA(cq, cr, i, w); cpairs[0] += 1
                        if a != b: cbad[0] += 1; print('CHILD DIM MISMATCH', sorted(cq.edges(keys=True)), i, w, a, b)
            if kd != kj: Jbad += 1; print('J MISMATCH', sorted(q.edges(keys=True)), R, v, kd, kj)
    print('quivers with a doubled arrow: %d; (i,v) pairs: %d; dim mismatches: %d; vertices tested for J: %d (J>0 in %d); J mismatches: %d\nmutated children (acyclic, gate-admitted, J=0): %d, (i,v) pairs %d, dim mismatches %d' % (npar, npairs, bad, nJ, Jpos, Jbad, nchild[0], cpairs[0], cbad[0]))

if MYARGS[1] == '--controls': run_controls()
elif MYARGS[1] == '--dimcheck': dimcheck(int(MYARGS[2]) if len(MYARGS) > 2 else 200, int(MYARGS[3]) if len(MYARGS) > 3 else 1)
