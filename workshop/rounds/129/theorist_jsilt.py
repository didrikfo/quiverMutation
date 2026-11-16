"""Round 129 (theorist), T10 mechanism: is every key-kept J != 0 step of E-323 a non-tilting irreducible silting mutation?

For a step (parent P, vertex k, side F|R) put a = P (F) or a = P^op (R) and build, in K^b(proj a), the Okuyama-Rickard complex
    T = T^+_k = Y + sum_{i != k} P_i,     Y = [ P_k --f--> sum_{beta: k -> t} P_t ]   (P_k in degree -1),
f = the arrows out of k.  Hom(P_x, P_y) = paths x ~> y mod I, composition = concatenation (Ladkani's p |-> p beta).  Computed
independently of `perI` / `isTiltingAtVertex`, by chain maps modulo homotopy (exact rationals):
    H[m]  = dim Hom(T, T[m]) for m = -2..2      (silting: H[m] = 0, m > 0; tilting: also H[m] = 0, m < 0)
    appr  = f is a left add(a/P_k)-approximation (every path k ~> j, j != k, factors through f) -- AI Def 2.30 input
    End0  = the matrix dim Hom(T_x, T_y) (Cartan of End(T) in degree 0), compared with the exact Cartan of the walk's child
            (child^op for R) up to simultaneous permutation.
Also the library's J (sum of `perI`-style kernel dims, own copy) and `mutation.isTiltingAtVertex`.
Controls (same test): E-080 (abde = acde, at d, F): J != 0 known; abd = acd (at d, F): tilting J = 0 known; the LNA J = 0 step of E-323's
crossing walk (first step, if any is J = 0); planted errors: Y built with one arrow dropped (not an approximation: H[1] must be > 0);
End0 compared with the child of a different vertex (must mismatch).
usage (repo root):  .venv/bin/python -W ignore workshop/rounds/129/theorist_jsilt.py [controls|steps|all]
"""
import sys, json, itertools, collections, warnings
from fractions import Fraction
sys.path.insert(0, '.')
import networkx as nx
from quivermutation import arrowPaths as ap, procedure, pathAlgebra, mutation, reduction
warnings.simplefilter('ignore')

# ---------- algebras ----------
def deser(d):
    A = pathAlgebra.PathAlgebra(); Q = nx.MultiDiGraph()
    for x, y, k in d['arrows']: Q.add_edge(x, y, key=k)
    A.quiver = Q
    ar = [ap.combination({tuple(tuple(a) for a in p): c for p, c in r}) for r in d['rels']]
    A.arrowRels = ar; A.rels = [[list(ap.projectPath(p)) for p in r] for r in ar]
    return A

def build(arrows, rels):
    """arrows: list (x, y); rels: list of {tuple of arrow indices: coef}"""
    arr = [[x, y, 0] for x, y in arrows]
    return deser(dict(arrows=arr, rels=[[[[arr[i] for i in p], c] for p, c in r.items()] for r in rels]))

class Cat:
    """the path category of a = kQ/I (acyclic): Hom(P_x, P_y) = paths x ~> y mod I."""
    def __init__(s, A):
        s.Q = A.quiver; s.R = procedure.relationsFrom(A); s._b = {}; s._p = {}
        assert nx.is_directed_acyclic_graph(s.Q), 'acyclic only'
    def piv(s, x, y):
        if (x, y) not in s._p: s._p[(x, y)] = ap.idealBasis(s.Q, s.R, x, y)
        return s._p[(x, y)]
    def basis(s, x, y):
        if (x, y) not in s._b:
            s._b[(x, y)] = [()] if x == y else [p for p in ap.allPathsBetween(s.Q, x, y) if p not in s.piv(x, y)]
        return s._b[(x, y)]
    def nf(s, comb, x, y):     # normal form as {basis path: coef}
        if x == y: return {(): Fraction(comb.get((), 0))} if comb.get((), 0) else {}
        return {k: Fraction(v) for k, v in ap.reduceAgainstPivots(ap.combination({p: c for p, c in comb.items()}), s.piv(x, y)).items()} if comb else {}
    def mul(s, f, g, x, y, z):  # f: x ~> y then g: y ~> z
        out = collections.Counter()
        for p, a in f.items():
            for q, b in g.items(): out[p + q] += a * b
        return s.nf({k: v for k, v in out.items() if v}, x, z)

def rank(rows):
    piv = {}; r = 0
    for row in rows:
        row = {k: Fraction(v) for k, v in row.items() if v}
        while row:
            h = min(row)
            if h not in piv: piv[h] = {k: v / row[h] for k, v in row.items()}; r += 1; break
            f, pv = row[h], piv[h]
            row = {k: row.get(k, 0) - f * pv.get(k, 0) for k in set(row) | set(pv)}; row = {k: v for k, v in row.items() if v}
    return r

# complex: dict deg -> list of vertices; diff: dict deg -> {(i, j): comb}  (summand i of deg -> summand j of deg + 1)
def homdim(C, cx, D, m):
    """dim Hom_K(X, Y[m]) for complexes X = cx[0], cx[1] and Y = D."""
    (X, dX), (Y, dY) = cx, D
    var = [(d, i, j, p) for d in X for i, x in enumerate(X[d]) for j, y in enumerate(Y.get(d + m, [])) for p in C.basis(x, y)]
    vix = {v: n for n, v in enumerate(var)}
    def compose_rows(gvec_of):   # unused helper placeholder
        pass
    # cycle condition: for each d: dX^d then g^{d+1}  ==  g^d then dY^{d+m}, as map X_d -> Y_{d+m+1}
    rows = collections.defaultdict(dict)      # condition key -> {var index: coef}
    for (d, i, j, p), n in vix.items():
        x, y = X[d][i], Y[d + m][j]
        for (j2, t), c in ((jt, c) for jt, c in dY.get(d + m, {}).items() if jt[0] == j):   # g then dY
            z = Y[d + m + 1][t]
            for q, a in C.mul({p: 1}, c, x, y, z).items(): rows[(d, i, t, q)][n] = rows[(d, i, t, q)].get(n, 0) + a
        for (i0, i2), c in dX.get(d - 1, {}).items():                                       # dX then g (g at degree d, source d-1)
            if i2 != i: continue
            w = X[d - 1][i0]
            for q, a in C.mul(c, {p: 1}, w, x, y).items(): rows[(d - 1, i0, j, q)][n] = rows[(d - 1, i0, j, q)].get(n, 0) - a
    # nullspace dimension of the condition matrix (conditions as linear forms in the variables)
    cond = [r for r in rows.values() if r]
    Z = len(var) - rank(cond)
    # boundaries: h^d: X_d -> Y_{d+m-1};  dh = dX^d then h^{d+1}  +  h^d then dY^{d+m-1}
    bnd = []
    for d in X:
        for i, x in enumerate(X[d]):
            for j, y in enumerate(Y.get(d + m - 1, [])):
                for p in C.basis(x, y):
                    vec = collections.Counter()
                    for (j0, t), c in dY.get(d + m - 1, {}).items():          # h^d then dY: X_d -> Y_{d+m}
                        if j0 != j: continue
                        for q, a in C.mul({p: 1}, c, x, y, Y[d + m][t]).items(): vec[vix[(d, i, t, q)]] += a
                    for (i0, i2), c in dX.get(d - 1, {}).items():             # dX^{d-1} then h^d: X_{d-1} -> Y_{d+m-1}
                        if i2 != i: continue
                        for q, a in C.mul(c, {p: 1}, X[d - 1][i0], x, y).items(): vec[vix[(d - 1, i0, j, q)]] += a
                    bnd.append(dict(vec))
    return Z - rank(bnd)

def T_plus(A, k, drop=0):
    """summands of T^+_k: list of (name, complex). drop = number of out-arrows omitted (planted error)."""
    outs = ap.arrowsOutOf(A.quiver, k)[drop:]
    Y = ({-1: [k], 0: [b[1] for b in outs]}, {-1: {(0, j): {(b,): 1} for j, b in enumerate(outs)}})
    return [('Y', Y)] + [(i, ({0: [i]}, {})) for i in sorted(A.quiver.nodes) if i != k]

def approx(C, A, k):
    """every element of Hom(P_k, P_j), j != k, factors through f (sum over arrows beta: k -> t of beta * Hom(P_t, P_j))."""
    outs = ap.arrowsOutOf(A.quiver, k)
    for j in A.quiver.nodes:
        if j == k: continue
        B = C.basis(k, j)
        if not B: continue
        img = [C.mul({(b,): 1}, {q: 1}, k, b[1], j) for b in outs for q in C.basis(b[1], j)]
        if rank(img) < len(B): return False
    return True

def J_lib(A, k):   # own copy of round 033 perI: sum of kernel dims of p |-> (p beta)_beta
    C = Cat(A); outs = ap.arrowsOutOf(A.quiver, k); tot = 0
    for i in A.quiver.nodes:
        if i == k: continue
        B = C.basis(i, k)
        if not B: continue
        rows = [{(b, q): v for b in outs for q, v in C.mul({p: 1}, {(b,): 1}, i, k, b[1]).items()} for p in B]
        tot += len(B) - rank(rows)
    return tot

def cartan(A):
    C = Cat(A); V = sorted(A.quiver.nodes)
    return [[len(C.basis(x, y)) for y in V] for x in V]

def perm_eq(M, N):
    def g(M):
        G = nx.DiGraph()
        for a in range(len(M)):
            G.add_node(a, w=M[a][a])
            for b in range(len(M)):
                if a != b and M[a][b]: G.add_edge(a, b, w=M[a][b])
        return G
    if len(M) != len(N) or sorted(sum(M, [])) != sorted(sum(N, [])): return False
    return nx.is_isomorphic(g(M), g(N), node_match=lambda u, v: u['w'] == v['w'], edge_match=lambda u, v: u['w'] == v['w'])

def test(A, k, child, drop=0):
    """A = the algebra mutated (P or P^op), child = the child on the same side (child or child^op)."""
    C = Cat(A); T = T_plus(A, k, drop)
    H = {m: sum(homdim(C, X[1], Y[1], m) for X in T for Y in T) for m in (-2, -1, 1, 2)}
    Hm1_YY = homdim(C, T[0][1], T[0][1], -1)
    E0 = [[homdim(C, X[1], Y[1], 0) for Y in T] for X in T]
    Cc = cartan(child) if child is not None else None
    return dict(J=J_lib(A, k), lib_tilt=mutation.isTiltingAtVertex(A, k), appr=approx(C, A, k), H=H, Hm1_YY=Hm1_YY,
                silting=all(H[m] == 0 for m in (1, 2)), tilting=all(H[m] == 0 for m in (-2, -1, 1, 2)),
                End0_eq_child=None if Cc is None else perm_eq(E0, Cc), End0_eq_childT=None if Cc is None else perm_eq(E0, [list(r) for r in zip(*Cc)]),
                endY=E0[0][0])

def show(name, r):
    print('%-34s J %d lib_tilt %-5s appr %-5s H[-2,-1,1,2] %s %s %s %s  silting %-5s tilting %-5s End0~child %s (~child^T %s) dimEnd(Y) %d' % (
        name, r['J'], r['lib_tilt'], r['appr'], r['H'][-2], r['H'][-1], r['H'][1], r['H'][2], r['silting'], r['tilting'],
        r['End0_eq_child'], r['End0_eq_childT'], r['endY']), flush=True)

def child_of(A, k):
    raw = mutation.quiverMutationAtVertex(A, k)
    if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): return None
    return reduction.reducePathAlgebra(raw)

def controls():
    sq = [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]          # a: 1->2, b: 1->3, c: 2->4, d': 3->4, e: 4->5 ; E-080 numbering 1..5
    E080 = build(sq, [{(0, 2, 4): 1, (1, 3, 4): -1}])     # path 1-2-4-5 = 1-3-4-5, mutate at 4
    SQ = build(sq, [{(0, 2): 1, (1, 3): -1}])              # 1-2-4 = 1-3-4, mutate at 4
    for name, A in (('E-080 abde=acde @4 F (J!=0 known)', E080), ('abd=acd @4 F (tilting known)', SQ)):
        show(name, test(A, 4, child_of(A, 4)))
    show('PLANTED drop arrow: 1->2->{3,4} @2', test(build([(1, 2), (2, 3), (2, 4)], []), 2, None, drop=1))
    show('PLANTED drop arrow: E-080 @1', test(E080, 1, None, drop=1))
    A = build([(1, 2), (2, 3), (3, 4), (4, 5)], [{(0, 1, 2): 1}])   # LNA A5 with 1-2-3-4 = 0
    for k in (2, 3):
        show('LNA A5/(abc) @%d F' % k, test(A, k, child_of(A, k)))
    show('PLANTED wrong child: LNA @2 vs child@3', test(A, 2, child_of(A, 3)))

def steps():
    rows = [json.loads(l) for l in open('workshop/rounds/125/toolsmith_jxclass_out.jsonl')]
    real = [r for r in rows if 'distinct' not in r and not r['planted']]
    tab = collections.Counter(); first = True
    for r in real:
        P, Ch = deser(r['parent_struct']), deser(r['child_struct']); k = int(r['v'][1:])
        if r['v'][0] == 'F': A, CA = P, Ch
        else: A, CA = pathAlgebra.dualPathAlgebra(P), pathAlgebra.dualPathAlgebra(Ch)
        t = test(A, k, CA)
        show('%s %s|%s %s' % (r['root'], r['side'], r['path'], r['v']), t)
        tab[(r['v'][0], t['J'] > 0, t['lib_tilt'], t['appr'], t['silting'], t['tilting'], t['H'][-1] == t['J'], t['End0_eq_child'], t['End0_eq_childT'], t['Hm1_YY'])] += 1
        if first and r['root'] == '0040033030':
            first = False
    print('(side, J>0, lib_tilt, approximation, silting, tilting, H[-1] == J, End0 ~ child, End0 ~ child^T, dim Hom(Y, Y[-1])) -> count')
    for k, v in sorted(tab.items(), key=str): print('  ', k, v)

def walk(root='0040033030', side='root', path='F3 R9 F3 F3'):
    """replay one recorded crossing walk from its LNA root and run the test on EVERY step (J = 0 ones are the in-family tilting control)."""
    from quivermutation import nakayama as nk
    X = nk.LinearNakayamaAlgebra.fromClassName(root)
    if side == 'dual': X = pathAlgebra.dualPathAlgebra(X)
    for st in path.split():
        k = int(st[1:]); A = X if st[0] == 'F' else pathAlgebra.dualPathAlgebra(X)
        assert mutation.mutationIsPossibleAtVertex(A, k)
        ch = child_of(A, k)
        show('walk %s %s: %s' % (root, side, st), test(A, k, ch))
        X = ch if st[0] == 'F' else pathAlgebra.dualPathAlgebra(ch)

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if mode in ('controls', 'all'): controls()
    if mode in ('controls', 'all'): walk(); walk('003033030', 'root', 'F3 R8 F3 F3')
    if mode in ('steps', 'all'): steps()
