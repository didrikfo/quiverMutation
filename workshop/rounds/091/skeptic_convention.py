"""Round 091 (skeptic, T5): left/right and op convention of the repo's F = mutateAtVertex and R = mutateLeftAtVertex
against INDEPENDENT silting/tilting mutation in K^b(proj A) of a monomial algebra A (zero relations).

Own engine (nothing from the package except the objects audited and the quiver/Cartan readers):
 * convention fixed once: P_v = projective REPRESENTATION at v (paths starting at v); Hom(P_u, P_w) has the nonzero paths w ~> u as basis;
   a morphism given by the path p is "precomposition"; composition P_u -> P_v -> P_w is the path (w~>v) followed by (v~>u).
   This is the convention of arXiv:2112.08129 Sect. 2 (P_i = paths starting in i, P_j -> P_i for a path i ~> j).
 * T_R = A/P_i + X_R, X_R = cocone of the right add(A/P_i)-approximation sum_{alpha out of i} P_{t(alpha)} -> P_i, i.e. the complex
   [sum P_{t alpha} (deg 0) -> P_i (deg 1)]  (AI mu^- , 2112.08129 "right tilting mutation").
 * T_L = A/P_i + X_L, X_L = cone of the left approximation P_i -> sum_{beta into i} P_{s(beta)}, i.e. [P_i (deg -1) -> sum P_{s beta} (deg 0)] (AI mu^+).
 * M_X[x][y] = dim Hom_K(T_x, T_y);  A_X[x][y] = dim of the part of Hom(T_x,T_y) not a sum of composites through the other summands (Gabriel arrows x -> y of End(T) read as maps T_x -> T_y).
For a child (quiver, relations) with N[a][b] = number of nonzero paths a ~> b (cartanMatrix stores matrix[target][source]) and arrows G[a][b]:
 at T = A the labelled identity is N[a][b] = dim Hom(P_b, P_a) = M[b][a], so "child = End(T) in this convention" means N = M^T, G = A^T.
usage: skeptic_convention.py demo | sweep N_VERTICES SEEDS
"""
import sys, random, itertools
sys.path.insert(0, '.')
P = 32003

# ---------- linear algebra mod P ----------
def rref(rows, ncols):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        k = next((k for k in range(r, len(rows)) if rows[k][c] % P), None)
        if k is None: continue
        rows[r], rows[k] = rows[k], rows[r]
        inv = pow(rows[r][c], P - 2, P)
        rows[r] = [x * inv % P for x in rows[r]]
        for k in range(len(rows)):
            if k != r and rows[k][c] % P:
                f = rows[k][c]; rows[k] = [(x - f * y) % P for x, y in zip(rows[k], rows[r])]
        piv.append(c); r += 1
        if r == len(rows): break
    return rows[:r], piv

def nullspace(M, ncols):           # vectors v with M v = 0 ; M list of rows
    if not M: return [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    R, piv = rref(M, ncols); free = [c for c in range(ncols) if c not in piv]; out = []
    for f in free:
        v = [0] * ncols; v[f] = 1
        for r, c in zip(R, piv): v[c] = (-r[f]) % P
        out.append(v)
    return out

def rank(rows, ncols): return len(rref(rows, ncols)[1]) if rows else 0

# ---------- the monomial algebra ----------
class Alg:
    def __init__(s, n, arrows, zero):          # arrows: list of (s,t); zero: list of vertex sequences (relations)
        s.n = n; s.arrows = [(a, b, k) for k, (a, b) in enumerate(arrows)]
        s.zero = [tuple(z) for z in zero]; s._paths = {}
    def vz(s, path):                           # path as tuple of arrow indices -> vertex sequence
        if not path: return ()
        a = s.arrows; return (a[path[0]][0],) + tuple(a[k][1] for k in path)
    def nonzero(s, path):
        seq = s.vz(path)
        for z in s.zero:
            for st in range(len(seq) - len(z) + 1):
                if seq[st:st + len(z)] == z: return False
        return True
    def paths(s, u, w):                        # nonzero paths u ~> w (tuples of arrow indices)
        if (u, w) in s._paths: return s._paths[(u, w)]
        out = []
        def rec(v, p):
            if v == w: out.append(tuple(p))
            for k, (a, b, _) in enumerate(s.arrows):
                if a == v and s.nonzero(tuple(p) + (k,)): rec(b, p + [k])
        rec(u, []); s._paths[(u, w)] = out; return out

# ---------- complexes of projectives ----------
class Cx:   # terms: {deg: [vertex,...]}; d: {deg: {(a,b): {path: coef}}} map term a (deg n) -> term b (deg n+1), path = w ~> u for P_u -> P_w
    def __init__(s, terms, d): s.terms = terms; s.d = d
def stalk(v): return Cx({0: [v]}, {})

def homspace(A, U, V, m):                      # basis of Hom^m(U,V)
    B = []
    for n, us in U.terms.items():
        for a, u in enumerate(us):
            for b, w in enumerate(V.terms.get(n + m, [])):
                for p in A.paths(w, u): B.append((n, a, b, p))
    return B

def compose_path(A, p_after, p_before):        # g after f: g has path w~>v, f has path v~>u -> w~>u
    q = p_after + p_before
    return q if A.nonzero(q) else None

def Dmat(A, U, V, m):
    """matrix of D: Hom^m -> Hom^{m+1}, D f = d_V f - (-1)^m f d_U, as a list of rows (dim target x dim source)"""
    S = homspace(A, U, V, m); T = homspace(A, U, V, m + 1)
    ti = {b: i for i, b in enumerate(T)}
    cols = []
    for (n, a, b, p) in S:
        col = {}
        # d_V after f : U^n -> V^{n+m} -> V^{n+m+1}
        for (b1, b2), comb in V.d.get(n + m, {}).items():
            if b1 != b: continue
            for q, c in comb.items():
                r = compose_path(A, q, p)
                if r is not None:
                    key = ti[(n, a, b2, r)]; col[key] = (col.get(key, 0) + c) % P
        # f after d_U : U^{n-1} -> U^n -> V^{n+m}   gives component at degree n-1
        for (a0, a1), comb in U.d.get(n - 1, {}).items():
            if a1 != a: continue
            for q, c in comb.items():
                r = compose_path(A, p, q)
                if r is not None:
                    key = ti[(n - 1, a0, b, r)]; col[key] = (col.get(key, 0) - (1 if m % 2 == 0 else -1) * c) % P
        cols.append(col)
    M = [[cols[j].get(i, 0) for j in range(len(S))] for i in range(len(T))]
    return M, S, T

class HomK:
    """H^m(Hom(U,V)) with explicit representative cochains and a coordinate map"""
    def __init__(s, A, U, V, m):
        s.A, s.U, s.V, s.m = A, U, V, m
        Dm, S, T = Dmat(A, U, V, m); s.S = S
        s.cyc = nullspace(Dm, len(S)) if S else []
        Dp, Sp, Tp = Dmat(A, U, V, m - 1)
        # boundaries: columns of Dp (images), as vectors in S
        s.bnd = [[Dp[i][j] for i in range(len(S))] for j in range(len(Sp))] if Sp else []
        s.bnd = [b for b in s.bnd if any(b)]
        Rb, pb = rref(s.bnd, len(S)) if s.bnd else ([], [])
        s.Rb, s.pb = Rb, pb
        # basis of quotient: extend boundary space by cycles
        cur = Rb[:]; s.basis = []
        for z in s.cyc:
            if rank(cur + [z], len(S)) > len(cur):
                cur = cur + [z]; s.basis.append(z)
        s.dim = len(s.basis)
        s.full = Rb + s.basis
        s.Fr, s.Fpiv = (rref(s.full, len(S)) if s.full else ([], []))
    def coords(s, v):
        """coordinates of the cycle v in the basis of the quotient"""
        nb = len(s.Rb); k = len(s.full)
        # solve sum x_i full_i = v
        if k == 0: return []
        cols = len(s.S)
        aug = [[s.full[i][j] for i in range(k)] + [v[j]] for j in range(cols)]
        R, piv = rref(aug, k + 1)
        assert k not in piv, 'not a cycle in span'
        x = [0] * k
        for r, c in zip(R, piv): x[c] = r[k]
        return x[nb:]

def cx_hom(A, U, V, m=0): return HomK(A, U, V, m)

def comp_classes(A, hUV, hVW, hUW):
    """matrix of composition Hom(U,V) x Hom(V,W) -> Hom(U,W): returns table[i][j] = coords of (g_j . f_i)"""
    tab = {}
    idx = {b: i for i, b in enumerate(hUW.S)}
    for i, f in enumerate(hUV.basis):
        for j, g in enumerate(hVW.basis):
            out = [0] * len(hUW.S)
            for ff, cf in zip(hUV.S, f):
                if not cf: continue
                for gg, cg in zip(hVW.S, g):
                    if not cg: continue
                    (n, a, b, p), (n2, b2, c, q) = ff, gg
                    if n2 != n or b2 != b: continue
                    r = compose_path(A, q, p)
                    if r is not None:
                        k = idx[(n, a, c, r)]; out[k] = (out[k] + cf * cg) % P
            tab[(i, j)] = hUW.coords(out)
    return tab

def mutate_cx(A, i, side):
    if side == 'R':   # cocone of sum_{alpha: i->j} P_j -> P_i : deg 0 : P_j's, deg 1 : P_i
        al = [(k, b) for k, (a, b, _) in enumerate(A.arrows) if a == i]
        X = Cx({0: [b for _, b in al], 1: [i]}, {0: {(a, 0): {(k,): 1} for a, (k, _) in enumerate(al)}})
    else:             # cone of P_i -> sum_{beta: h->i} P_h : deg -1 : P_i, deg 0 : P_h's  (map P_i -> P_h = path h ~> i = (beta,))
        be = [(k, a) for k, (a, b, _) in enumerate(A.arrows) if b == i]
        X = Cx({-1: [i], 0: [a for _, a in be]}, {-1: {(0, a): {(k,): 1} for a, (k, _) in enumerate(be)}})
    return X

def invariants(A, i, side, shifts=(-2, -1, 1, 2)):
    summ = [(v, stalk(v)) for v in range(A.n) if v != i] + [(i, mutate_cx(A, i, side))]
    summ.sort(key=lambda t: t[0])
    n = A.n
    H = {(x, y): cx_hom(A, summ[x][1], summ[y][1], 0) for x in range(n) for y in range(n)}
    vanish = True
    for m in shifts:
        for x in range(n):
            for y in range(n):
                if cx_hom(A, summ[x][1], summ[y][1], m).dim: vanish = False
    M = [[H[(x, y)].dim for y in range(n)] for x in range(n)]
    # arrows: Hom(x,y) modulo composites through other summands c != x,y
    Ar = [[0] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if x == y: continue
            imgs = []
            for c in range(n):
                if c in (x, y): continue
                tab = comp_classes(A, H[(x, c)], H[(c, y)], H[(x, y)])
                imgs += list(tab.values())
            Ar[x][y] = M[x][y] - rank(imgs, M[x][y]) if M[x][y] else 0
    return vanish, M, Ar

# ---------- the repo side ----------
def repo_children(A, i):
    from quivermutation import procedure, arrowPaths as ap
    from tests.helpers import path_algebra, quiet
    pa = path_algebra(arrows=[(a + 1, b + 1) for a, b, _ in A.arrows], rels=[[tuple(v + 1 for v in z)] for z in A.zero])
    q = pa.quiver; rels = procedure.relationsFrom(pa)
    res = {}
    for name, f in (('F', procedure.mutateAtVertex), ('R', procedure.mutateLeftAtVertex)):
        try:
            nq, nr = quiet(f, q, rels, i + 1)
            nq, nr = quiet(procedure.reduce, nq, nr)
        except Exception as e:
            res[name] = ('err', repr(e)[:80]); continue
        C = ap.cartanMatrix(nq, nr, exact=True)
        n = A.n
        N = [[C[b][a] for b in range(n)] for a in range(n)]            # N[a][b] = paths a ~> b
        G = [[0] * n for _ in range(n)]
        for t, h, k in nq.edges(keys=True): G[t - 1][h - 1] += 1
        res[name] = (N, G, ap.hasOrientedCycle(nq))
    return res, procedure.isMutable(q, rels, i + 1)

def T(M): return [list(r) for r in zip(*M)]

def verdict(child, inv):
    """which of the 4 candidates (M, M^T) x ... the child matches; child=(N,G), inv=(vanish,M,Ar)"""
    N, G, cyc = child; _, M, Ar = inv
    gm = (G == T(Ar)); nm = None if cyc else (N == T(M))
    return {'N=M^T,G=A^T': gm and (nm is not False), 'G=A^T': gm, 'N=M^T': nm, 'cyc': cyc,
            'N=M,G=A': (N == M and G == Ar) if not cyc else (G == Ar)}

def opposite(A):
    return Alg(A.n, [(b, a) for a, b, _ in A.arrows], [tuple(reversed(z)) for z in A.zero])

def opp_literal(A, i, kids):
    """Oppermann read literally (his P_i -> sum_{alpha out of i} P_{t alpha} has Hom(P_i,P_j) = paths i ~> j, i.e. A^op in my covariant engine):
    his LEFT silting mutation = my 'L' on A^op ; his RIGHT = my 'R' on A^op.  The child must equal End (no op): N = M', G = A'."""
    Ao = opposite(A); out = {}
    for name, side in (('F', 'L'), ('R', 'R')):
        if kids[name][0] == 'err': continue
        vanish, M, Ar = invariants(Ao, i, side)
        N, G, cyc = kids[name]
        out[name] = (vanish, G == Ar and (cyc or N == M), G == Ar and (cyc or N == T(M)))
    return out

def run_case(A, i):
    kids, gate = repo_children(A, i)
    out = {'gate': gate}
    invR = invariants(A, i, 'R'); invL = invariants(A, i, 'L')
    out['tiltR'], out['tiltL'] = invR[0], invL[0]
    out['R_vs_L_differ'] = (invR[1], invR[2]) != (invL[1], invL[2])
    for name in 'FR':
        if kids[name][0] == 'err': out[name] = 'err'; continue
        out[name] = {'vs T_R': verdict(kids[name], invR), 'vs T_L': verdict(kids[name], invL)}
    out['opp'] = opp_literal(A, i, kids)
    return out, kids, invR, invL

def demo():
    # non-symmetric: 0->1->2, 1->3, 0->4 ; zero relation 0->1->2.   vertices 0..4
    A = Alg(5, [(0, 1), (1, 2), (1, 3), (0, 4)], [(0, 1, 2)])
    for i in (0, 1, 2, 3, 4):
        out, kids, invR, invL = run_case(A, i)
        print('i =', i, out)
        if i == 0:
            print(' M_R', invR[1], '\n A_R', invR[2], '\n M_L', invL[1], '\n child F N', kids['F'][0], '\n child R N', kids['R'][0])

def random_alg(n, rng):
    arrows = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < 0.38]
    if len(arrows) < n - 1 or set(v for a in arrows for v in a) != set(range(n)): return None
    A0 = Alg(n, arrows, [])
    cand = []
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if (a, b) in arrows and (b, c) in arrows: cand.append((a, b, c))
    zero = [z for z in cand if rng.random() < 0.4]
    return Alg(n, arrows, zero)

def sweep(n, seeds):
    from collections import Counter
    cnt = Counter(); exs = []
    for seed in range(seeds):
        rng = random.Random(seed); A = random_alg(n, rng)
        if A is None: continue
        for i in range(n):
            if not any(a == i for a, _, _ in A.arrows): continue
            try: out, kids, invR, invL = run_case(A, i)
            except Exception as e:
                cnt['crash ' + repr(e)[:50]] += 1; exs.append(('crash', seed, i)); continue
            cnt['F cases: gate %s, T_R tilting %s' % (out['gate'], out['tiltR'])] += 1
            hasin = any(b == i for _, b, _ in A.arrows)
            for name, tilt, key in (('F', out['tiltR'], 'vs T_R'), ('R', out['tiltL'] and hasin, 'vs T_L')):
                if out[name] == 'err' or not tilt: cnt[name + ' not tilting/err'] += 1; continue
                other = 'vs T_L' if key == 'vs T_R' else 'vs T_R'
                v = out[name]
                cnt['%s tilting: child cyclic=%s ; own T: G=A^T %s ; N=M^T %s' % (name, v[key]['cyc'], v[key]['G=A^T'], v[key]['N=M^T'])] += 1
                cnt['%s tilting: other reading (child=End(T), no transpose) matches: %s' % (name, v[key]['N=M,G=A'])] += 1
                ov = out['opp'][name]
                cnt['%s Oppermann-literal (A^op, his %s): tilting %s ; child = End (no op) %s ; child = End^T reading %s' % (name, 'left' if name == 'F' else 'right', ov[0], ov[1], ov[2])] += 1
                if out['R_vs_L_differ']:
                    cnt['%s power case (T_R,T_L data differ): matches OTHER T (N=M^T,G=A^T): %s' % (name, v[other]['N=M^T,G=A^T'])] += 1
                    cnt['%s power case: own %s' % (name, v[key]['N=M^T,G=A^T'])] += 1
                    if not v[key]['N=M^T,G=A^T']: exs.append((seed, i, name))
    for k in sorted(cnt): print(cnt[k], k)
    print('failures to match own:', exs[:10])

if __name__ == '__main__':
    if sys.argv[1] == 'demo': demo()
    else: sweep(int(sys.argv[2]), int(sys.argv[3]))
