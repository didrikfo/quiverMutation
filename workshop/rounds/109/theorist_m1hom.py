"""Round 109 theorist, T5/H-023: Hom-computed End(T) for the algebra A = kQ/(a then b) of model M1 (E-284):
Q: a: i->j, b: j->i, relation 'a b' (a then b) = 0.  A has dim 5 (e_i, e_j, a, b, ba), gl.dim 2, Nakayama Kupisch (2,3), Ext^2(S_i,S_i)=1
(this is the degree-1 loop f of M1).  We form the two silting mutations of A at i (left: P_i -> P_j via b; right: P_j -> P_i via a),
check Hom(T,T[m]) for all m, and compute End(T) (dim, Cartan matrix, quiver by radical layers) and the Coxeter polynomial.
Right modules, P_x = e_x A = paths starting at x, Hom(P_x,P_y) = paths y->x, composition = concatenation.  Exact over GF(p) (p = 1000003).
Controls: (C1) a planted wrong approximation (P_i -> P_j via 0) must fail silting; (C2) the non-mutation T = A[1]+ P_j ... see code.
Usage (repo root): .venv/bin/python workshop/rounds/109/theorist_m1hom.py"""
import itertools, sys
P = 1000003
V = {'i': 0, 'j': 1}
# paths: name -> (src, tgt)
ARROWS = {'a': (0, 1), 'b': (1, 0)}
def paths(maxlen=6):
    ps = {(): None}
    out = []
    for v in (0, 1): out.append(((v,), v, v))     # trivial path e_v encoded as ('e',v)
    return out
# Represent a path as a tuple of arrow names; trivial path at v as ('e', v).  Relation: 'a','b' consecutive = 0 (a then b).
def src(p): return p[1] if p[0] == 'e' else ARROWS[p[0]][0]
def tgt(p): return p[1] if p[0] == 'e' else ARROWS[p[-1]][1]
def norm(p):
    if p[0] == 'e': return p
    for k in range(len(p) - 1):
        if p[k] == 'a' and p[k + 1] == 'b': return None
    return p
def cat(p, q):      # p then q (concatenation); None if zero or incomposable
    if tgt(p) != src(q): return None
    if p[0] == 'e': return q
    if q[0] == 'e': return p
    return norm(p + q)
ALL = [('e', 0), ('e', 1)]
frontier = [(x,) for x in ARROWS]
while frontier:
    nxt = []
    for p in frontier:
        if norm(p) is None: continue
        ALL.append(p)
        for x in ARROWS:
            if ARROWS[p[-1]][1] == ARROWS[x][0]: nxt.append(p + (x,))
    frontier = nxt
assert len(ALL) == 5, ALL
def hom_basis(x, y):   # Hom(P_x, P_y) = paths y -> x
    return [p for p in ALL if src(p) == y and tgt(p) == x]

# A complex of projectives: dict degree -> list of vertices (summands); differential d[k]: matrix entries (row summand of X^{k+1}, col summand of X^k) = dict path->coeff, as a map X^k -> X^{k+1}
class Cx:
    def __init__(self, terms, diffs): self.terms = terms; self.diffs = diffs   # diffs[k] = {(r,c): {path: coeff}}
def shift(X, s):   # X[s]: (X[s])^k = X^{k+s}
    return Cx({k - s: v for k, v in X.terms.items()}, {k - s: d for k, d in X.diffs.items()})
def dsum(X, Y):
    terms = {}; diffs = {}
    for k in set(X.terms) | set(Y.terms):
        terms[k] = [('X', v) for v in X.terms.get(k, [])] + [('Y', v) for v in Y.terms.get(k, [])]
    for k in set(X.diffs) | set(Y.diffs):
        dd = {}
        nx = len(X.terms.get(k, [])); ny1 = len(X.terms.get(k + 1, []))
        for (r, c), e in X.diffs.get(k, {}).items(): dd[(r, c)] = e
        for (r, c), e in Y.diffs.get(k, {}).items(): dd[(r + ny1, c + nx)] = e
        diffs[k] = dd
    cx = Cx({k: [v for _, v in t] for k, t in terms.items()}, diffs)
    return cx
def mult(e1, e2):  # e1 then e2 with element dicts (the composite of maps: first e2 ... careful)  composition g∘f for f: P_x->P_y via z (y->x), g: P_y->P_w via z' (w->y): z' then z
    out = {}
    for p, c in e1.items():
        for q, d in e2.items():
            r = cat(p, q)
            if r is not None: out[r] = (out.get(r, 0) + c * d) % P
    return {k: v for k, v in out.items() if v}

def hom_complex(X, Y, m):
    """Hom^m(X,Y) = prod_k Hom(X^k, Y^{k+m}); returns basis list of (k, r, c, path) with r summand of Y^{k+m}, c summand of X^k"""
    B = []
    for k, xs in X.terms.items():
        ys = Y.terms.get(k + m, [])
        for c, xv in enumerate(xs):
            for r, yv in enumerate(ys):
                for p in hom_basis(xv, yv): B.append((k, r, c, p))
    return B
def apply_D(X, Y, m, f):     # f: dict basis->coeff in Hom^m ; D f = d_Y f - (-1)^m f d_X  in Hom^{m+1}
    out = {}
    def add(key, c): out[key] = (out.get(key, 0) + c) % P
    for (k, r, c, p), co in f.items():
        # d_Y f : X^k -> Y^{k+m} -> Y^{k+m+1}: entries (r2, r) of Y.diffs[k+m] composed with f(r,c): composite = (path of d) after (path of f): f path z (y->x), d path z' (w->y)  => z' then z
        for (r2, r1), e in Y.diffs.get(k + m, {}).items():
            if r1 != r: continue
            for q, d in e.items():
                t = cat(q, p)
                if t is not None: add((k, r2, c, t), co * d)
        # f d_X : X^{k-1} -> X^k -> Y^{k+m}: d_X entry (c, c0) path q (x_c -> x_c0 ... as path from X^k summand c to X^{k-1} summand c0), then f path p (y -> x_c): composite p then q
        sgn = -1 if m % 2 == 0 else 1
        for (c1, c0), e in X.diffs.get(k - 1, {}).items():
            if c1 != c: continue
            for q, d in e.items():
                t = cat(p, q)
                if t is not None: add((k - 1, r, c0, t), sgn * co * d)
    return {k: v for k, v in out.items() if v}

def rank_mod(rows, ncols):
    rows = [list(r) for r in rows]; rk = 0; col = 0
    piv_rows = []
    for col in range(ncols):
        pr = None
        for i in range(rk, len(rows)):
            if rows[i][col] % P: pr = i; break
        if pr is None: continue
        rows[rk], rows[pr] = rows[pr], rows[rk]
        inv = pow(rows[rk][col], P - 2, P)
        rows[rk] = [x * inv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col] % P:
                f = rows[i][col]; rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def Dmatrix(X, Y, m):
    src_b = hom_complex(X, Y, m); tgt_b = hom_complex(X, Y, m + 1); idx = {b: n for n, b in enumerate(tgt_b)}
    rows = []
    for b in src_b:
        img = apply_D(X, Y, m, {b: 1}); row = [0] * len(tgt_b)
        for key, c in img.items(): row[idx[key]] = c
        rows.append(row)
    return rows, len(src_b), len(tgt_b)
def hom_dim(X, Y, m):
    """dim Hom_{K^b(proj)}(X, Y[m]) = dim H^m Hom(X,Y)"""
    rows_m, n, nt = Dmatrix(X, Y, m); rk_m = rank_mod(rows_m, nt) if rows_m and nt else 0
    rows_p, n2, nt2 = Dmatrix(X, Y, m - 1); rk_p = rank_mod(rows_p, nt2) if rows_p and nt2 else 0
    return n - rk_m - rk_p
def d2check(X):
    for k in X.diffs:
        for (r2, r1), e in X.diffs.get(k + 1, {}).items():
            pass
    return True

def proj(v, k=0): return Cx({k: [v]}, {})
def cone_map(src_v, tgt_v, path, deg_src):
    """two-term complex  P_src (degree deg_src) --path--> P_tgt (degree deg_src+1)"""
    return Cx({deg_src: [src_v], deg_src + 1: [tgt_v]}, {deg_src: {(0, 0): {path: 1}}})

def analyse(name, T, Tsum_names, verbose=True):
    # T is a list of complexes (summands); returns Hom table
    rng = range(-3, 4)
    tab = {}
    for m in rng:
        tab[m] = [[hom_dim(T[x], T[y], m) for y in range(len(T))] for x in range(len(T))]
    nonzero = {m: tab[m] for m in rng if m != 0 and any(any(r) for r in tab[m])}
    print(name, ': summands', Tsum_names)
    print('  Hom(T_x, T_y) at m = 0 (row x, col y):', tab[0])
    print('  m != 0 with Hom(T,T[m]) != 0:', {m: v for m, v in nonzero.items()} if nonzero else 'none (tilting, in the range -3..3)')
    return tab

if __name__ == '__main__':
    Pj, Pi = proj(1), proj(0)
    A = [Pi, Pj]
    print('A itself:'); tabA = analyse('A', A, ['P_i', 'P_j'])
    # left mutation: X = (P_i --b--> P_j), P_i in degree 0 (so that X is the cone of P_i -> P_j placed with P_j in degree 1 ... both placements tried)
    for label, path, s, t_, deg in [('left  P_i --b--> P_j', ('b',), 0, 1, 0), ('right P_j --a--> P_i', ('a',), 1, 0, 0)]:
        for sh in (-1, 0, 1):
            X = shift(cone_map(s, t_, path, deg), sh)
            T = [X, Pj] if s == 0 else [X, Pj]
            print('\n', label, 'shift', sh)
            tab = analyse('T', T, ['X', 'P_j'])
    print('\nC1 planted error: map P_i -> P_j by the zero map (X = P_i[1] + P_j): not an approximation')
    X = Cx({0: [0, 1]}, {})   # P_i + P_j in degree 0 (zero differential) -- cone of the zero map P_i->P_j with deg shift ignored
    X0 = Cx({0: [0], 1: [1]}, {0: {}})   # zero differential P_i -> P_j
    analyse('C1', [X0, Pj], ['X0', 'P_j'])
