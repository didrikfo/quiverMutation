"""Round 107 (toolsmith), H-027 / T10: breadth-first search of the single-summand EXCHANGE graph of tilting complexes of K^b(proj B_n),
starting at T0 = B_n (= +P_i), through TILTING complexes only (every node passes Hom(T,T[m]) = 0, m != 0, by the round-050 chain-map Hom,
and has |V| summands with K0 classes independent).  Nodes are sets of (r, a) = X_r[a], X_r the n(n-1) indecomposables found by
toolsmith_cones.py.  Hom tables h[r][s][j] = dim Hom(X_r, X_s[j]) computed once by homdim.  Edge: two nodes sharing |V|-1 summands.
Per node: Cartan H(End T)[x][y] = dim Hom(Z_x, Z_y); test H ~ Cartan(A_n) up to simultaneous permutation / transpose.
Controls: (i) depth 1 reproduces E-282's 8 (n=6) / 6 (n=5) tilting one-step complexes among the vertex-local ones, with the 2 silting-only
ones rejected; (ii) planted: an object set with a Hom(T,T[1]) entry is rejected.
Usage: toolsmith_tiltsearch.py n depth [window]"""
import sys, itertools, json, time
sys.path.insert(0, '.')
exec(open('workshop/rounds/107/toolsmith_cones.py').read().split("if __name__ == '__main__'")[0])
import numpy as np
from quivermutation import invariants
from collections import Counter

def build(n):
    alg = Bn(n); A = Alg(alg)
    objs = [deser(S) for S in json.load(open('workshop/rounds/107/toolsmith_objs_%d.json' % n))]
    N = len(objs); J = 14
    h = [[{j: homdim(A, objs[r], shiftC(objs[s], j), 0) for j in range(-J, J + 1)} for s in range(N)] for r in range(N)]
    return A, objs, h, J

def perm_canon_eq(X, Y):
    n = len(X)
    if sorted(X.flatten()) != sorted(Y.flatten()): return False
    rs = lambda M: sorted(map(tuple, M.tolist()))
    for p in itertools.permutations(range(n)):
        if (X[np.ix_(p, p)] == Y).all(): return True
    return False

if __name__ == '__main__':
    n = int(sys.argv[1]); depth = int(sys.argv[2]); W = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    t0 = time.time()
    A, objs, h, J = build(n); N = len(objs)
    print('n', n, 'objects up to shift', N, 'Hom tables %.0fs' % (time.time() - t0)); sys.stdout.flush()
    cls = [np.array(klassC(X, n)) for X in objs]
    def hom(x, y, m=0):          # dim Hom(Z_x, Z_y[m]), Z=(r,a)=X_r[a]
        (r, a), (s, b) = x, y; j = b - a + m
        return h[r][s][j] if -J <= j <= J else 0
    def compat(x, y): return all(hom(x, y, m) == 0 and hom(y, x, m) == 0 for m in range(-2 * W - 2, 2 * W + 3) if m != 0)
    cache = {}
    def cp(x, y):
        if (x, y) not in cache: cache[(x, y)] = compat(x, y)
        return cache[(x, y)]
    def selfok(x): return all(hom(x, x, m) == 0 for m in range(-8, 9) if m != 0) and hom(x, x, 0) == 1
    def klass(x): return ((-1) ** x[1]) * cls[x[0]]
    # projectives: base object whose layers are a single vertex i at degree 0
    P = {}
    for r, X in enumerate(objs):
        if len(X[0]) == 1 and len(X[0][0]) == 1: P[X[0][0][0]] = (r, 0)
    T0 = frozenset(P[i] for i in range(1, n + 1))
    assert len(T0) == n, P
    allz = [(r, a) for r in range(N) for a in range(-W, W + 1)]
    def isTilting(T):
        T = list(T)
        return all(cp(x, y) for x, y in itertools.combinations(T, 2)) and abs(round(np.linalg.det(np.array([klass(x) for x in T])))) == 1
    assert isTilting(T0)
    def H(T): T = sorted(T); return np.array([[hom(x, y) for y in T] for x in T])
    CA = np.array(invariants.cartanMatrix(An(n), exact=True).tolist(), dtype=int)
    sigA = (sorted(CA.sum(0)), sorted(CA.sum(1)), sorted(CA.diagonal()))
    def isA(T):
        M = H(T)
        if sorted(M.flatten()) != sorted(CA.flatten()): return False
        return perm_canon_eq(M, CA) or perm_canon_eq(M.T, CA)
    def neighbours(T):
        out = []
        for U in T:
            rest = [x for x in T if x != U]
            for y in allz:
                if y == U or y in rest or not selfok(y): continue
                if all(cp(y, x) for x in rest):
                    T2 = frozenset(rest + [y])
                    if abs(round(np.linalg.det(np.array([klass(x) for x in T2])))) == 1: out.append((U, y, T2))
        return out
    seen = {T0: 0}; frontier = [T0]; parent = {T0: None}; hitsA = []
    print('depth 0 : 1 node; End T0 Cartan == A_n?', isA(T0))
    for dpt in range(1, depth + 1):
        nxt = []
        for T in frontier:
            for U, y, T2 in neighbours(T):
                if T2 not in seen:
                    seen[T2] = dpt; parent[T2] = (T, U, y); nxt.append(T2)
        hit = [T for T in nxt if isA(T)]
        print('depth %d : %d new tilting complexes (cumulative %d); Cartan(End T) ~ Cartan(A_%d): %d   [%.0fs]' % (dpt, len(nxt), len(seen), n, len(hit), time.time() - t0)); sys.stdout.flush()
        hitsA += [(dpt, T) for T in hit]; frontier = nxt
        if not nxt: break
    def show(T):
        return ' + '.join('X%d[%d]' % x for x in sorted(T))
    print('first hits:'); 
    for dpt, T in hitsA[:4]:
        path = []; c = T
        while parent[c]: path.append(parent[c][1:]); c = parent[c][0]
        print(' depth', dpt, show(T), ' path (replaced -> new):', list(reversed(path))); print('   H =', H(T).tolist())
    # control: one-step neighbours of T0 with their H, and the gate-style (vertex-local) cut
    print('== one-step neighbours of B_%d (all, not only vertex-local), tilting by construction ==' % n)
    for U, y, T2 in neighbours(T0):
        print(' replace', U, '->', y, 'class', klass(y).tolist(), 'Cartan~A_n', isA(T2))
    # silting-only controls: for each (summand U of T0, base object y) that is exceptional but NOT compatible with the other n-1 summands
    bad = 0; tot = 0
    for U in T0:
        rest = [x for x in T0 if x != U]
        for y in allz:
            if y in T0 or not selfok(y): continue
            ok = all(cp(y, x) for x in rest); tot += 1
            if not ok: bad += 1
    print('planted/non-tilting control: %d of %d (exceptional y, removed summand) pairs from B_%d are rejected as non-tilting' % (bad, tot, n))
    json.dump({'n': n, 'depth': depth, 'nodes': len(seen), 'hits': [[d, sorted(T)] for d, T in hitsA]}, open('workshop/rounds/107/toolsmith_tiltsearch_%d.json' % n, 'w'))

# ---------------------------------------------------------------------------------------------------------------
# Strong check on hits: End(T) ~ A_n as ALGEBRA, not only its Cartan matrix.  All Hom(Z_x,Z_y) are <= 1-dim (Cartan 0/1).  A basic algebra with
# Cartan matrix = Cartan(A_n) up to perm is A_n iff the composite of the n-2 consecutive "arrow" maps is non-zero in homotopy classes
# (then every sub-composite is non-zero by associativity; the (1,n) entry is 0).  Done on explicit chain maps.
def homotopyImage(A, X, Y):
    """list of chain maps d s + s d for basis homotopies s: X^k -> Y^{k-1}."""
    Xd, Xdi = X; Yd, Ydi = Y; out = []
    for k in sorted(set(Xd)):
        if k - 1 not in Yd: continue
        for a, va in enumerate(Xd[k]):
            for b, vb in enumerate(Yd[k - 1]):
                for p in A.basis(va, vb):
                    g = {}
                    # d_Y s : X^k -> Y^k
                    for (b1, c), cm in Ydi.get(k - 1, {}).items():
                        if b1 != b: continue
                        vc = Yd[k][c]
                        comb = {p + pc: Fraction(cc) for pc, cc in cm.items()}
                        for pp, cf in A.nf(va, vc, comb).items(): g.setdefault(k, {}).setdefault((a, c), {})[pp] = g.get(k, {}).get((a, c), {}).get(pp, 0) + cf
                    # s d_X : X^{k-1} -> Y^{k-1}
                    for (a0, a1), cm in Xdi.get(k - 1, {}).items():
                        if a1 != a: continue
                        v0 = Xd[k - 1][a0]
                        comb = {pc + p: Fraction(cc) for pc, cc in cm.items()}
                        for pp, cf in A.nf(v0, vb, comb).items(): g.setdefault(k - 1, {}).setdefault((a0, b), {})[pp] = g.get(k - 1, {}).get((a0, b), {}).get(pp, 0) + cf
                    out.append(g)
    return out
def flat(g): return {(k, ab, p): v for k, m in g.items() for ab, c in m.items() for p, v in c.items() if v}
def rk(vecs):
    keys = {}; rows = []
    for v in vecs:
        rows.append({keys.setdefault(k, len(keys)): Fraction(x) for k, x in v.items()})
    piv = {}; r = 0
    for row in rows:
        row = dict(row)
        while row:
            h = min(row)
            if h in piv:
                f = row[h]
                for u, x in piv[h].items(): row[u] = row.get(u, 0) - f * x
                row = {u: x for u, x in row.items() if x}
            else:
                inv = 1 / row[h]; piv[h] = {u: x * inv for u, x in row.items()}; r += 1; break
    return r
def generator(A, X, Y):
    """a chain map X->Y that is not null-homotopic (None if Hom = 0)."""
    cm = chainMaps(A, X, Y); hi = [flat(g) for g in homotopyImage(A, X, Y)]; base = rk(hi)
    for g in cm:
        if rk(hi + [flat(g)]) > base: return g
    return None
def compose(A, X, Y, Z, f, g):   # g o f
    out = {}
    for k in f:
        if k not in g: continue
        for (a, b), fc in f[k].items():
            for (b2, c), gc in g[k].items():
                if b2 != b: continue
                comb = out.setdefault(k, {}).setdefault((a, c), {})
                for p1, v1 in fc.items():
                    for p2, v2 in gc.items(): comb[p1 + p2] = comb.get(p1 + p2, 0) + v1 * v2
    res = {}
    for k, m in out.items():
        for (a, c), comb in m.items():
            nf = A.nf(X[0][k][a], Z[0][k][c], {q: Fraction(v) for q, v in comb.items()})
            if nf: res.setdefault(k, {})[(a, c)] = nf
    return res
def isNullHomotopic(A, X, Y, g):
    hi = [flat(x) for x in homotopyImage(A, X, Y)]
    return rk(hi + [flat(g)]) == rk(hi)
def verifyAlgebraA(A, T, n):
    """returns (True, order) if End(T) ~ A_n as an algebra (given Cartan ~ A_n), else (False, reason)."""
    T = sorted(T); Zs = [shiftC(objs[r], a) for r, a in T]
    M = H(T); CAm = CA
    for transp in (False, True):
        MM = M.T if transp else M
        for p in itertools.permutations(range(n)):
            if (MM[np.ix_(p, p)] == CAm).all():
                order = [T[i] for i in p]; Zo = [Zs[i] for i in p]
                # in CAm, nonzero (i,j) for i<=j: Hom(Z_i, Z_j) (or transposed: Hom(Z_j, Z_i)); chain of arrows i -> i+1
                chain = []
                for i in range(n - 1):
                    X, Y = (Zo[i + 1], Zo[i]) if not transp else (Zo[i], Zo[i + 1])
                    chain.append((X, Y, generator(A, X, Y)))
                if any(g is None for _, _, g in chain): return False, 'missing arrow map'
                res = []
                for s in (0, 1):   # composites of consecutive maps 0..n-3 and 1..n-2 (length n-2)
                    seq = list(reversed(chain[s:s + n - 2])) if not transp else chain[s:s + n - 2]
                    X0 = seq[0][0]; cur = seq[0][2]; Y0 = seq[0][1]
                    for X, Y, g in seq[1:]:
                        cur = compose(A, X0, Y0, Y, cur, g); Y0 = Y
                    res.append(bool(cur) and not isNullHomotopic(A, X0, Y0, cur))
                return all(res), ('perm %s transposed %s composites nonzero %s' % (p, transp, res))
    return False, 'no perm'
print('== algebra check on hits: End(T) ~ A_n as algebra (composites of consecutive arrow maps nonzero in homotopy classes) ==')
okh = 0
for dpt, T in hitsA[:int(sys.argv[4]) if len(sys.argv) > 4 else 6]:
    ok, msg = verifyAlgebraA(A, T, n)
    okh += ok; print(' depth', dpt, show(T), '->', ok, msg)
mx = max(abs(a) for T in seen for r, a in T); print('max |shift| among', len(seen), 'nodes:', mx, '(window', W, ')')
print('== controls for the algebra check ==')
def projC(i): return ({0: [i]}, {})
for nm, alg_, vs in (('A_%d, n-2 arrows (must be NONZERO)' % n, An(n), list(range(1, n))), ('A_%d, n-1 arrows (must VANISH)' % n, An(n), list(range(1, n + 1))), ('B_%d, 2 arrows (must be NONZERO)' % n, Bn(n), list(range(n - 3, n))), ('B_%d (3-arrow relation: composite of 3 consecutive arrows must vanish)' % n, Bn(n), list(range(n - 3, n + 1)))):
    AA = Alg(alg_); done = False
    for direction in (0, 1):
        ch = []
        for u, v in zip(vs, vs[1:]):
            X, Y = (projC(u), projC(v)) if direction == 0 else (projC(v), projC(u))
            ch.append((X, Y, generator(AA, X, Y)))
        if any(g is None for _, _, g in ch): continue
        seq = ch if direction == 0 else list(reversed(ch))
        X0, Y0, cur = seq[0]
        for X, Y, g in seq[1:]:
            cur = compose(AA, X0, Y0, Y, cur, g); Y0 = Y
        nz = bool(cur) and not isNullHomotopic(AA, X0, Y0, cur)
        print(' ', nm, ': composite of', len(vs) - 1, 'consecutive arrows nonzero =', nz); done = True; break
    if not done: print(' ', nm, ': arrow maps not found')
