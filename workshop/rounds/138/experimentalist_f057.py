"""Round 138 (experimentalist): independent rebuild of F-057 (a), the n = 11 pure-F crossing.

Recorded path (E-332, rounds/129/maverick.md): dual(303305000) . F6 F11 F7 F11 F7 F7* F3 F3 F4 F7 -> LNA 403033030
(* = the key-kept J != 0 step at depth 6); second root 303500030 . F1 F1 F4 F3 F1 F8*.

Only the MUTATION itself is taken from the library (mutation.quiverMutationAtVertex, reduction.reducePathAlgebra,
mutationIsPossibleAtVertex, dualPathAlgebra, arrowPaths.isIllegalRelation): that is the object under test.
Everything that certifies is written here: LNA decoding and Cartan matrix from the name, the Cartan matrix of a quiver with
relations by exact linear algebra over Q (paths modulo the two-sided ideal span, no monomial assumption), Coxeter
polynomial (Faddeev-LeVerrier over Fractions), Smith normal form (own elimination, checked against sympy), the n = 11
key-group census from an own enumeration of LNAs, and the identification of the end node as an LNA.
No round-122/125/127/129 script is imported or exec'd.  The library's integerCartanMatrix / coxeter key are called ONLY
as a cross-check column (`lib=`), never in a verdict.

usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/138/experimentalist_f057.py [walk|census|controls|all]
"""
import sys, itertools, random, time
from fractions import Fraction
sys.path.insert(0, '.')

N = 11

# ---------------------------------------------------------------- own linear algebra
def rank_q(rows, ncols):
    """rank over Q of a list of dict{col: Fraction}"""
    piv = {}  # col -> row(dict) with leading col
    r = 0
    for row in rows:
        row = {c: Fraction(v) for c, v in row.items() if v}
        while row:
            c = min(row)
            if c in piv:
                p = piv[c]; f = row[c] / p[c]
                for cc, vv in p.items():
                    nv = row.get(cc, 0) - f * vv
                    if nv: row[cc] = nv
                    else: row.pop(cc, None)
            else:
                piv[c] = row; r += 1; break
    return r

def rank_p(rows, p):
    """rank over GF(p) of a list of dict{col: int}"""
    piv = {}; r = 0
    for row in rows:
        row = {c: v % p for c, v in row.items() if v % p}
        while row:
            c = min(row)
            if c in piv:
                q = piv[c]; f = row[c] * pow(q[c], -1, p) % p
                for cc, vv in q.items():
                    nv = (row.get(cc, 0) - f * vv) % p
                    if nv: row[cc] = nv
                    else: row.pop(cc, None)
            else:
                piv[c] = row; r += 1; break
    return r

def inv_q(M):
    n = len(M); A = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

def det_q(M):
    n = len(M); A = [[Fraction(x) for x in row] for row in M]; d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return 0
        if p != c: A[c], A[p] = A[p], A[c]; d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f: A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return d

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def charpoly(M):
    """Faddeev-LeVerrier; coefficients low degree first, monic"""
    n = len(M); M = [[Fraction(x) for x in r] for r in M]
    c = [Fraction(0)] * (n + 1); c[n] = Fraction(1)
    Mk = [[Fraction(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = [[Mk[i][j] + (c[n - k + 1] if i == j else 0) for j in range(n)] for i in range(n)]
        AM = matmul(M, Mk)
        c[n - k] = -sum(AM[i][i] for i in range(n)) / k
        Mk = AM
    assert all(x.denominator == 1 for x in c)
    return tuple(int(x) for x in c)

def coxeter_key(C):
    """charpoly of Phi = -C^T C^-1 (det C must be 1)"""
    assert det_q(C) == 1
    n = len(C); Ci = inv_q(C); CT = [[C[j][i] for j in range(n)] for i in range(n)]
    Phi = [[-x for x in row] for row in matmul(CT, Ci)]
    return charpoly(Phi)

def snf(M):
    """own Smith normal form: sorted absolute diagonal"""
    A = [list(map(int, r)) for r in M]; m, n = len(A), len(A[0]); diag = []
    t = 0
    while t < min(m, n):
        nz = [(abs(A[i][j]), i, j) for i in range(t, m) for j in range(t, n) if A[i][j]]
        if not nz: break
        _, i, j = min(nz)
        A[t], A[i] = A[i], A[t]
        for row in A: row[t], row[j] = row[j], row[t]
        while True:
            done = True
            for i in range(t + 1, m):
                if A[i][t]:
                    q = A[i][t] // A[t][t]; A[i] = [x - q * y for x, y in zip(A[i], A[t])]
                    if A[i][t]: done = False
            for j in range(t + 1, n):
                if A[t][j]:
                    q = A[t][j] // A[t][t]
                    for row in A: row[j] -= q * row[t]
                    if A[t][j]: done = False
            if done:
                bad = [(i, j) for i in range(t + 1, m) for j in range(t + 1, n) if A[i][j] % A[t][t]]
                if not bad: break
                i, j = bad[0]; A[t] = [x + y for x, y in zip(A[t], A[i])]; done = False
            if not done:
                nz = [(abs(A[i][j]), i, j) for i in range(t, m) for j in range(t, n) if A[i][j] and (i == t or j == t)]
                _, i, j = min(nz); A[t], A[i] = A[i], A[t]
                for row in A: row[t], row[j] = row[j], row[t]
        diag.append(abs(A[t][t])); t += 1
    diag += [0] * (min(m, n) - len(diag))
    return tuple(sorted(diag))

def minor_gcd(S, k):
    """gcd of all k x k minors (k = n - 1 here): a congruence invariant computed without any SNF"""
    from math import gcd
    n = len(S); g = 0
    for rs in itertools.combinations(range(n), k):
        for cs in itertools.combinations(range(n), k):
            g = gcd(g, int(det_q([[S[i][j] for j in cs] for i in rs])))
    return g

def sym(C):
    n = len(C); return [[C[i][j] + C[j][i] for j in range(n)] for i in range(n)]

def signature(S):
    import numpy as np
    ev = np.linalg.eigvalsh(np.array(S, dtype=float)); return int((ev > 1e-9).sum()), int((ev < -1e-9).sum())

# ---------------------------------------------------------------- own LNA code
def lna_intervals(name):
    n = len(name) + 2
    return n, sorted((i + 1, i + 1 + int(d)) for i, d in enumerate(name) if d != '0')

def lna_cartan(n, iv):
    """C[i][j] = dim e_i A e_j for the linear quiver 1 -> 2 -> ... -> n: 1 iff i <= j and no relation interval inside [i, j]"""
    return [[int(i <= j and not any(i <= a and b <= j for a, b in iv)) for j in range(1, n + 1)] for i in range(1, n + 1)]

def lna_name(n, iv):
    s = ['0'] * (n - 2)
    for a, b in iv: assert b - a <= 9; s[a - 1] = str(b - a)
    return ''.join(s)

def dual_iv(n, iv):
    return sorted((n + 1 - b, n + 1 - a) for a, b in iv)

def all_lnas(n):
    """antichains (no containment) of intervals [a, b], b - a >= 2: every LNA with no length-1 relations"""
    out = []
    def rec(a, last_b, cur):
        if a > n - 2: out.append(list(cur)); return
        rec(a + 1, last_b, cur)
        for b in range(max(a + 2, last_b + 1), n + 1):   # starts increase, so ends must strictly increase
            cur.append((a, b)); rec(a + 1, b, cur); cur.pop()
    rec(1, 0, []); return out

# ---------------------------------------------------------------- own Cartan of a quiver with relations
def all_paths(Q):
    """every path (tuple of edges (u, v, k)) incl. trivial ones as ('e', v); Q acyclic"""
    edges = sorted(Q.edges(keys=True)); out_e = {}
    for e in edges: out_e.setdefault(e[0], []).append(e)
    paths = [(v,) for v in sorted(Q.nodes)]
    stack = [(e,) for e in edges]
    while stack:
        p = stack.pop(); paths.append(p)
        for e in out_e.get(p[-1][1], []): stack.append(p + (e,))
    return paths

def src(p): return p[0] if not isinstance(p[0], tuple) else p[0][0]
def tgt(p): return p[0] if not isinstance(p[0], tuple) else p[-1][1]
def cat(p, q):
    a = () if not isinstance(p[0], tuple) else p
    b = () if not isinstance(q[0], tuple) else q
    return a + b if a + b else p

def cartan_general(A, relations, char=None):
    Q = A.quiver; import networkx as nx
    assert nx.is_directed_acyclic_graph(Q)
    V = sorted(Q.nodes); paths = all_paths(Q)
    idx = {p: i for i, p in enumerate(paths)}
    into = {v: [p for p in paths if tgt(p) == v] for v in V}
    outof = {v: [p for p in paths if src(p) == v] for v in V}
    gens = {}
    for r in relations:
        terms = [(tuple(p), c) for p, c in r.items() if c]
        s = {src(p) for p, _ in terms}; t = {tgt(p) for p, _ in terms}
        assert len(s) == 1 and len(t) == 1, 'inhomogeneous relation'
        s, t = s.pop(), t.pop()
        for p in into[s]:
            for q in outof[t]:
                row = {}
                for w, c in terms:
                    k = idx[cat(cat(p, w), q)]; row[k] = row.get(k, 0) + c
                gens.setdefault((s, t if True else None, src(p), tgt(q)), []).append(row)
    C = []
    for i in V:
        rowC = []
        for j in V:
            ps = [p for p in paths if src(p) == i and tgt(p) == j]
            rows = [g for (s, t, a, b), gl in gens.items() if a == i and b == j for g in gl]
            rowC.append(len(ps) - (rank_q(rows, len(paths)) if char is None else rank_p(rows, char)))
        C.append(rowC)
    return V, C

# ---------------------------------------------------------------- the walk (library mutation, own certificates)
def lib():
    from quivermutation import nakayama as nk, pathAlgebra, mutation, reduction, procedure, arrowPaths, invariants, search
    return nk, pathAlgebra, mutation, reduction, procedure, arrowPaths, invariants, search

def step(Alg, v):
    nk, pathAlgebra, mutation, reduction, procedure, ap, invariants, search = lib()
    assert mutation.mutationIsPossibleAtVertex(Alg, v), 'gate refuses F%d' % v
    raw = mutation.quiverMutationAtVertex(Alg, v)
    assert not any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)), 'illegal relation at F%d' % v
    return reduction.reducePathAlgebra(raw)

def certs(Alg):
    nk, pathAlgebra, mutation, reduction, procedure, ap, invariants, search = lib()
    rels = procedure.relationsFrom(Alg)
    V, C = cartan_general(Alg, rels)
    S = sym(C); d = det_q(C)
    charp = all(cartan_general(Alg, rels, q)[1] == C for q in (2, 3, 5))
    key = coxeter_key(C) if d == 1 else None
    libC = invariants.integerCartanMatrix(Alg)
    mono = all(sum(1 for c in r.values() if c) == 1 for r in rels)
    return dict(C=C, det=d, key=key, snf=snf(S), sig=signature(S), libagree=('C' if [list(map(int, r)) for r in libC] == C else
                     'C^T' if [list(map(int, r)) for r in zip(*libC)] == C else 'NO'),
                libkey=search._coxeterKeyOrNone(Alg), charp=charp,
                nrels=len(rels), mono=mono, arrows=Q_arrows(Alg))

def Q_arrows(Alg): return Alg.quiver.number_of_edges()

def as_lna(Alg):
    """if Alg is a linearly oriented A_n with monomial relations, return its LNA name (own reading); else None"""
    nk, pathAlgebra, mutation, reduction, procedure, ap, invariants, search = lib()
    Q = Alg.quiver; n = Q.number_of_nodes()
    if Q.number_of_edges() != n - 1 or any(Q.out_degree(v) > 1 or Q.in_degree(v) > 1 for v in Q.nodes): return None
    start = [v for v in Q.nodes if Q.in_degree(v) == 0]
    if len(start) != 1: return None
    order = [start[0]]
    while Q.out_degree(order[-1]): order.append(next(iter(Q.successors(order[-1]))))
    if len(order) != n: return None
    pos = {v: i + 1 for i, v in enumerate(order)}
    iv = []
    for r in procedure.relationsFrom(Alg):
        t = [p for p, c in r.items() if c]
        if len(t) != 1: return None
        p = t[0]; iv.append((pos[p[0][0]], pos[p[-1][1]]))
    iv = sorted(iv)
    # minimal generators only (a reduced presentation should already be minimal)
    iv = [x for x in iv if not any(y != x and x[0] <= y[0] and y[1] <= x[1] for y in iv)]
    return lna_name(n, iv), iv, order

def permuted(C, V, order):
    """C indexed by sorted labels V -> C indexed along the A_n path order"""
    ix = [V.index(v) for v in order]; return [[C[i][j] for j in ix] for i in ix]

def walk(root, side, moves, label):
    nk, pathAlgebra, mutation, reduction, procedure, ap, invariants, search = lib()
    r = nk.LinearNakayamaAlgebra.fromClassName(root)
    A = r if side == 'root' else pathAlgebra.dualPathAlgebra(r)
    # own check that the library's start object is the LNA we think it is
    n, iv = lna_intervals(root); iv0 = iv if side == 'root' else dual_iv(n, iv)
    got = as_lna(A); assert got and got[1] == iv0, ('start object mismatch', got, iv0)
    c0 = certs(A); assert permuted(c0['C'], sorted(A.quiver.nodes), got[2]) == lna_cartan(n, iv0)
    print('WALK %s: start %s(%s) = LNA %s, own Cartan == own LNA Cartan; SNF(S) %s sig %s key %s' %
          (label, side, root, lna_name(n, iv0), c0['snf'], c0['sig'], c0['key']))
    rows = [('0', 'start', c0)]; prev = c0
    for d, mv in enumerate(moves, 1):
        A = step(A, int(mv[1:])); c = certs(A)
        rows.append((str(d), mv, c))
        print('  d=%-2d %-4s arrows %2d rels %2d mono %-5s det %s keyKept %-5s SNF tail %s sig %s  lib Cartan %s, lib key == own %s, Cartan same over GF(2,3,5) %s%s' %
              (d, mv, c['arrows'], c['nrels'], c['mono'], c['det'], c['key'] == c0['key'], c['snf'][-3:], c['sig'], c['libagree'], c['libkey'] == c['key'], c['charp'],
               '   <-- SNF(S) CHANGES' if c['snf'] != prev['snf'] else ''), flush=True)
        prev = c
    end = as_lna(A)
    if end:
        en, eiv, eorder = end; Ce = lna_cartan(n, eiv)
        print('  end node reads as LNA %s (dual reading %s); own Cartan (path order) == LNA Cartan from that name: %s' %
              (en, lna_name(n, dual_iv(n, eiv)), Ce == permuted(prev['C'], sorted(A.quiver.nodes), eorder)))
        Ln = lna_cartan(n, eiv); kL = coxeter_key(Ln)
        print('  from scratch: LNA %s key == start key %s; SNF(S) %s vs start %s -> separated %s' %
              (en, kL == c0['key'], snf(sym(Ln))[-3:], c0['snf'][-3:], snf(sym(Ln)) != c0['snf']))
        S0 = sym(lna_cartan(n, iv0)); S1 = sym(Ln)
        print('  SNF-free check: det S = %s / %s; D_10 = gcd of all 10x10 minors of S = %s / %s (start / end)' %
              (det_q(S0), det_q(S1), minor_gcd(S0, n - 1), minor_gcd(S1, n - 1)))
    else:
        print('  end node is not a linear monomial A_n presentation')
    return rows, A

# ---------------------------------------------------------------- census of the key group
def fast_key(C):
    """integer Faddeev-LeVerrier in numpy int64 for a unitriangular upper C (the census); overflow-guarded.
    Group members are re-keyed with the exact Fraction coxeter_key below."""
    import numpy as np
    n = len(C); A = np.array(C, dtype=np.int64); Ci = np.eye(n, dtype=np.int64)
    for i in range(n - 1, -1, -1):          # back-substitution: A upper unitriangular
        for j in range(i + 1, n):
            if A[i, j]: Ci[i] -= A[i, j] * Ci[j]
    M = -A.T @ Ci; c = [0] * (n + 1); c[n] = 1; Mk = np.zeros((n, n), dtype=np.int64)
    for k in range(1, n + 1):
        Mk = Mk + c[n - k + 1] * np.eye(n, dtype=np.int64); AM = M @ Mk
        tr = int(np.trace(AM)); assert tr % k == 0 and np.abs(AM).max() < 2 ** 40
        c[n - k] = -tr // k; Mk = AM
    return tuple(c)

def census(rootkey_names):
    t0 = time.time(); L = all_lnas(N); groups = {}
    for iv in L:
        C = lna_cartan(N, iv); k = fast_key(C)
        groups.setdefault(k, []).append(iv)
    print('CENSUS n = %d: %d LNAs (no length-1 relations), %d keys, %.0fs' % (N, len(L), len(groups), time.time() - t0))
    out = {}
    for nm in rootkey_names:
        n, iv = lna_intervals(nm); k = coxeter_key(lna_cartan(n, iv)); G = groups[k]
        exact = sum(coxeter_key(lna_cartan(N, g)) == k for g in G)
        near = [g for g in L if g not in G][:2000:7]
        leak = sum(coxeter_key(lna_cartan(N, g)) == k for g in near)
        print('  exact Fraction key re-check: %d of %d group members equal; %d of %d sampled non-members equal (must be 0)'
              % (exact, len(G), leak, len(near)))
        names = {lna_name(N, g) for g in G}
        print('  duality-closed: %s; %d duality orbits' % (all(lna_name(N, dual_iv(N, g)) in names for g in G),
              len({min(lna_name(N, g), lna_name(N, dual_iv(N, g))) for g in G})))
        by = {}
        for g in G: by.setdefault(snf(sym(lna_cartan(N, g))), []).append(lna_name(N, g))
        print('  key of %s: group of %d LNAs; split by SNF(C+C^T): %s' % (nm, len(G), {s[-3:]: len(v) for s, v in by.items()}))
        out[nm] = (k, by)
        for s, v in sorted(by.items()):
            print('    SNF tail %s (%d): %s' % (s[-3:], len(v), ' '.join(sorted(v))))
    return groups, out

# ---------------------------------------------------------------- controls
def controls():
    import sympy
    from sympy.matrices.normalforms import smith_normal_form
    rnd = random.Random(138); bad = 0
    for t in range(300):
        m = rnd.randint(2, 7); M = [[rnd.randint(-6, 6) for _ in range(m)] for _ in range(m)]
        if t % 3 == 0: M[rnd.randrange(m)] = [0] * m
        D = smith_normal_form(sympy.Matrix(M), domain=sympy.ZZ)
        ref = tuple(sorted(abs(int(D[i, i])) for i in range(m)))
        if ref != snf(M): bad += 1; print('  SNF mismatch', M, ref, snf(M))
    print('CONTROL snf vs sympy: %d of 300 random integer matrices disagree' % bad)
    # congruence invariance (truth YES) and planted change (truth NO)
    n, iv = lna_intervals('403033030'); S = sym(lna_cartan(n, iv)); yes = 0; no = 0
    for t in range(50):
        U = [[int(i == j) for j in range(n)] for i in range(n)]
        for _ in range(30):
            i, j = rnd.sample(range(n), 2); c = rnd.choice([-1, 1])
            for r in range(n): U[r][j] += c * U[r][i]
        T = matmul(matmul([[U[j][i] for j in range(n)] for i in range(n)], S), U)
        yes += snf(T) == snf(S)
    print('CONTROL congruence U^T S U (truth: same SNF): %d of 50 same' % yes)
    # planted end-node error: the 13-class LNA 303500030 put in place of the end node must NOT separate from the root
    for a, b in [('303500030', '500033030'), ('003033030', '303500030'), ('403033030', '400330330'), ('003033030', '030330030')]:
        na, ia = lna_intervals(a); nb, ib = lna_intervals(b)
        Ca, Cb = lna_cartan(na, ia), lna_cartan(nb, ib)
        print('CONTROL pair %s / %s: key equal %s, SNF tails %s / %s -> separated %s' %
              (a, b, coxeter_key(Ca) == coxeter_key(Cb), snf(sym(Ca))[-3:], snf(sym(Cb))[-3:], snf(sym(Ca)) != snf(sym(Cb))))
    # planted Cartan error at the end node: change one entry -> key equality must fail (det kept = 1: upper entry)
    Cp = [row[:] for row in lna_cartan(n, iv)]; Cp[0][5] += 1
    print('CONTROL planted Cartan entry (end LNA, C[1][6] += 1): key still equal to start key %s' %
          (coxeter_key(Cp) == coxeter_key(lna_cartan(*lna_intervals('303305000')))))

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    t0 = time.time()
    if mode in ('walk', 'all'):
        walk('303305000', 'dual', 'F6 F11 F7 F11 F7 F7 F3 F3 F4 F7'.split(), 'b (join, 10 steps)')
        walk('303500030', 'root', 'F1 F1 F4 F3 F1 F8'.split(), 'a (crossing only, 6 steps)')
    if mode in ('census', 'all'):
        census(['303305000'])
    if mode in ('controls', 'all'):
        controls()
    print('total %.0fs' % (time.time() - t0))
