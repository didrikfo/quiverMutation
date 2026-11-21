"""Round 134, theorist, T5 / H-023: child = End(T) as ALGEBRAS for gate-admitted J = 0 steps, n <= 6 acyclic parents.

For a parent A = kQ/I (acyclic) and a gate-admitted vertex k, T = Y + sum_{i != k} P_i with Y = [P_k -> sum_{k->t} P_t]
(the Okuyama-Rickard complex; Hom engine of rounds/129/theorist_jsilt.py, independent of `perI` / `isTiltingAtVertex`).
  1. multiplication table of E = End_K(T) in degree 0: a basis of Z/B for every pair of summands, composition of chain
     maps, reduction mod boundaries (exact rationals).  Slot (x, y) = Hom(T_x, T_y), product f*g = "f then g"; labels
     Y <-> k, P_i <-> i.
  2. the library child C = reducePathAlgebra(quiverMutationAtVertex(A, k)), as quiver + relation generators; its
     dimensions from the round-131 engine (rounds/131/theorist_thin.py `Alg`, own row reduction, no library maths).
  3. ISO TEST (label-preserving, both orientations tried: C -> E and C^op -> E):
       (a) cheap invariants of the tables: Cartan, dim (rad/rad^2)_{xy} = number of child arrows x -> y;
       (b) exact: unknown images phi(a) = sum c_{a,b} b over a basis of rad(E)_{xy} for every child arrow a: x -> y;
           equations phi(r) = 0 for every relation generator r; Rabinowitsch 1 - t * prod det (images independent mod
           rad^2 in every slot); vertex-torus normalisation (one leading coefficient = 1 per edge of a spanning forest of
           the single-arrow slots: conjugation by sum t_x e_x is an automorphism of C).  Groebner basis over QQ (sympy):
           ISO <=> GB != [1].  ISO + equal dimension + generation mod rad^2 => phi is an algebra isomorphism C ~ E.
  Planted wrong children (same quiver, same Cartan matrix, different algebra), built from the true child:
     KILL  : a binomial relation p - l q replaced by the monomial p (and q kept nonzero): dims unchanged when the
             slot dimension is unchanged (checked), multiplication changed;
     SCALE : l -> 2 l (usually ISO, the torus absorbs it; the control that the test is not over-strict);
     OP    : the child read backwards is tried as the second orientation anyway.
usage (repo root):  .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py controls
                    .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py lna NMAX
                    .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py rand N SEED NMAX
"""
import sys, collections, random, itertools, time, signal, warnings
from fractions import Fraction
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/129'); sys.path.insert(0, 'workshop/rounds/131')
warnings.simplefilter('ignore')
import networkx as nx
import sympy as sp
import theorist_jsilt as J          # Cat, homdim, T_plus, J_lib, child_of, build
import theorist_thin as TH           # Alg (independent kQ/I engine), randomAlg, toLib
from quivermutation import mutation, procedure, pathAlgebra, arrowPaths as ap

# ---------------- exact linear algebra ----------------
def rref(rows, ncols):
    """rows: list of lists of Fractions. returns (pivot rows dict col->row, list of pivot cols)."""
    M = [list(r) for r in rows]; piv = []; r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; f = M[r][c]; M[r] = [x / f for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                g = M[i][c]; M[i] = [a - g * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    return M[:r], piv

def nullspace(rows, ncols):
    R, piv = rref(rows, ncols)
    free = [c for c in range(ncols) if c not in piv]; out = []
    for f in free:
        v = [Fraction(0)] * ncols; v[f] = Fraction(1)
        for row, c in zip(R, piv): v[c] = -row[f]
        out.append(v)
    return out

def solve_in(basis, target):
    """coefficients x with sum x_i basis_i = target (basis independent), or None."""
    n = len(basis); m = len(target)
    rows = [[basis[i][j] for i in range(n)] + [target[j]] for j in range(m)]
    R, piv = rref(rows, n + 1)
    if n in piv: return None
    x = [Fraction(0)] * n
    for row, c in zip(R, piv): x[c] = row[n]
    return x

# ---------------- End(T) multiplication table ----------------
class EndT:
    def __init__(s, A, k):
        s.C = J.Cat(A); s.k = k
        s.T = J.T_plus(A, k)                                    # [('Y', cx), (i, cx), ...]
        s.lab = [k if nm == 'Y' else nm for nm, _ in s.T]
        s.idx = {l: n for n, l in enumerate(s.lab)}
        s.V = {}; s.basis = {}; s.B = {}
        for X in range(len(s.T)):
            for Yy in range(len(s.T)): s._hom(X, Yy)
    def _vars(s, X, Yy):
        (Xc, _), (Yc, _) = s.T[X][1], s.T[Yy][1]
        return [(d, i, j, p) for d in Xc for i, x in enumerate(Xc[d]) for j, y in enumerate(Yc.get(d, [])) for p in s.C.basis(x, y)]
    def _hom(s, X, Yy):
        C = s.C; (Xc, dX), (Yc, dY) = s.T[X][1], s.T[Yy][1]
        var = s._vars(X, Yy); vix = {v: n for n, v in enumerate(var)}; N = len(var)
        rows = collections.defaultdict(lambda: [Fraction(0)] * N)
        for (d, i, j, p), n in vix.items():
            x, y = Xc[d][i], Yc[d][j]
            for (j2, t), c in dY.get(d, {}).items():
                if j2 != j: continue
                for q, a in C.mul({p: 1}, c, x, y, Yc[d + 1][t]).items(): rows[(d, i, t, q)][n] += a
            for (i0, i2), c in dX.get(d - 1, {}).items():
                if i2 != i: continue
                for q, a in C.mul(c, {p: 1}, Xc[d - 1][i0], x, y).items(): rows[(d - 1, i0, j, q)][n] -= a
        Z = nullspace(list(rows.values()), N) if rows else [[Fraction(int(a == b)) for a in range(N)] for b in range(N)]
        bnd = []
        for d in Xc:
            for i, x in enumerate(Xc[d]):
                for j, y in enumerate(Yc.get(d - 1, [])):
                    for p in C.basis(x, y):
                        v = [Fraction(0)] * N
                        for (j0, t), c in dY.get(d - 1, {}).items():
                            if j0 != j: continue
                            for q, a in C.mul({p: 1}, c, x, y, Yc[d][t]).items(): v[vix[(d, i, t, q)]] += a
                        for (i0, i2), c in dX.get(d - 1, {}).items():
                            if i2 != i: continue
                            for q, a in C.mul(c, {p: 1}, Xc[d - 1][i0], x, y).items(): v[vix[(d - 1, i0, j, q)]] += a
                        bnd.append(v)
        Bb, _ = rref(bnd, N) if bnd else ([], [])
        # complement of B inside Z
        basis = []; cur = list(Bb)
        for z in Z:
            if len(rref(cur + [z], N)[1]) > len(cur): basis.append(z); cur = rref(cur + [z], N)[0]
        s.V[(X, Yy)] = (var, vix); s.basis[(X, Yy)] = basis; s.B[(X, Yy)] = Bb
    def coords(s, X, Z, vec):
        x = solve_in(s.basis[(X, Z)] + s.B[(X, Z)], vec)
        assert x is not None, 'composite not a cycle?'
        return x[:len(s.basis[(X, Z)])]
    def compose(s, X, Yy, Z, f, g):
        """f in Hom(X,Y), g in Hom(Y,Z) as vectors over var(X,Y), var(Y,Z); returns vector over var(X,Z)."""
        C = s.C; (Xc, _), (Yc, _), (Zc, _) = s.T[X][1], s.T[Yy][1], s.T[Z][1]
        vXY, _ = s.V[(X, Yy)]; vYZ, _ = s.V[(Yy, Z)]; _, iXZ = s.V[(X, Z)]
        out = [Fraction(0)] * len(iXZ)
        gi = collections.defaultdict(list)
        for n, (d, j, l, q) in enumerate(vYZ):
            if g[n]: gi[(d, j)].append((l, q, g[n]))
        for n, (d, i, j, p) in enumerate(vXY):
            if not f[n]: continue
            for l, q, b in gi[(d, j)]:
                for r, c in C.mul({p: 1}, {q: 1}, Xc[d][i], Yc[d][j], Zc[d][l]).items():
                    out[iXZ[(d, i, l, r)]] += f[n] * b * c
        return out
    def table(s):
        """structure constants M[(X,Y,Z)][a][b] = coords of basis_a(X,Y) * basis_b(Y,Z) in basis(X,Z)."""
        s.M = {}
        n = len(s.T)
        for X, Yy, Z in itertools.product(range(n), repeat=3):
            s.M[(X, Yy, Z)] = [[s.coords(X, Z, s.compose(X, Yy, Z, f, g)) for g in s.basis[(Yy, Z)]] for f in s.basis[(X, Yy)]]
        # unit of each e_X: the identity chain map
        s.one = {}
        for X in range(n):
            var, vix = s.V[(X, X)]
            v = [Fraction(int(p == () and i == j)) for (d, i, j, p) in var]
            s.one[X] = s.coords(X, X, v)
        return s

    # radical and rad^2 per slot (in basis coordinates)
    def rad(s):
        n = len(s.T); s.R = {}
        for X in range(n):
            for Yy in range(n):
                d = len(s.basis[(X, Yy)])
                if X != Yy: s.R[(X, Yy)] = [[Fraction(int(a == b)) for a in range(d)] for b in range(d)]
                else:   # local ring: rad = kernel of trace of left multiplication
                    tr = [sum(s.M[(X, X, X)][a][b][b] for b in range(d)) for a in range(d)]
                    s.R[(X, Yy)] = nullspace([tr], d)
        s.R2 = {}
        for X in range(n):
            for Z in range(n):
                vecs = []
                for Yy in range(n):
                    for u in s.R[(X, Yy)]:
                        for w in s.R[(Yy, Z)]:
                            vecs.append(s.mulvec(X, Yy, Z, u, w))
                d = len(s.basis[(X, Z)])
                s.R2[(X, Z)] = rref(vecs, d)[0] if vecs else []
        return s
    def mulvec(s, X, Yy, Z, u, w):
        d = len(s.basis[(X, Z)]); out = [Fraction(0)] * d
        M = s.M[(X, Yy, Z)]
        for a, ua in enumerate(u):
            if not ua: continue
            for b, wb in enumerate(w):
                if not wb: continue
                for c in range(d): out[c] += ua * wb * M[a][b][c]
        return out
    def cartan(s):
        return {(s.lab[X], s.lab[Yy]): len(s.basis[(X, Yy)]) for X in range(len(s.T)) for Yy in range(len(s.T))}
    def irr(s):
        return {(s.lab[X], s.lab[Yy]): len(s.R[(X, Yy)]) - len(s.R2[(X, Yy)]) for X in range(len(s.T)) for Yy in range(len(s.T))}

# ---------------- the child as quiver + relations ----------------
def child_data(ch):
    arrows = [(x, y, kk) for x, y, kk in ch.quiver.edges(keys=True)]
    rels = [dict(r) for r in procedure.relationsFrom(ch)] if False else [dict(r) for r in ch.arrowRels]
    nodes = sorted(ch.quiver.nodes)
    return nodes, arrows, rels

def alg_cartan(nodes, arrows, rels):
    al = TH.Alg(len(nodes), arrows, rels)
    return {(x, y): (1 if x == y else len(al.std(x, y))) for x in nodes for y in nodes}

def reverse(nodes, arrows, rels):
    ra = [(y, x, kk) for x, y, kk in arrows]; m = {a: (a[1], a[0], a[2]) for a in arrows}
    return nodes, ra, [{tuple(m[a] for a in reversed(p)): c for p, c in r.items()} for r in rels]

class Timeout(Exception): pass
def _alarm(sig, frm): raise Timeout()

def iso(E, nodes, arrows, rels, tmax=90):
    """label-preserving: does there exist an algebra iso kQ_c/(rels) -> E?  returns 'ISO', 'NO:<reason>', 'UNDECIDED'."""
    lab = E.lab; ix = E.idx
    if sorted(lab) != sorted(nodes): return 'NO:vertices'
    if E.cartan() != alg_cartan(nodes, arrows, rels): return 'NO:cartan'
    cnt = collections.Counter((a[0], a[1]) for a in arrows)
    irr = E.irr()
    if any(irr[(x, y)] != cnt.get((x, y), 0) for x in nodes for y in nodes): return 'NO:irr'
    # unknowns
    syms = {}; img = {}
    for a in arrows:
        X, Yy = ix[a[0]], ix[a[1]]
        R = E.R[(X, Yy)]
        cs = [sp.Symbol('c%d_%d' % (arrows.index(a), b)) for b in range(len(R))]
        syms[a] = cs
        d = len(E.basis[(X, Yy)])
        img[a] = [sum((cs[b] * R[b][c] for b in range(len(R))), sp.Integer(0)) for c in range(d)]
    # leading coordinates: complement of rad^2 in rad, per slot
    dets = []; subs = {}
    slots = collections.defaultdict(list)
    for a in arrows: slots[(a[0], a[1])].append(a)
    forest = nx.Graph(); forest.add_nodes_from(nodes)
    for (x, y), al in slots.items():
        X, Yy = ix[x], ix[y]; d = len(E.basis[(X, Yy)])
        R2 = E.R2[(X, Yy)]
        # linear functionals vanishing on rad^2: nullspace of R2 (as rows), restricted to rad
        fun = nullspace(R2, d) if R2 else [[Fraction(int(a == b)) for a in range(d)] for b in range(d)]
        # keep the functionals that are independent on rad
        Rr = E.R[(X, Yy)]
        vals = [[sum(f[c] * r[c] for c in range(d)) for r in Rr] for f in fun]
        Rv, pv = rref(vals, len(Rr)); fun2 = Rv[:len(al)]
        assert len(fun2) == len(al)
        Mx = sp.Matrix([[sum((syms[a][b] * fr[b] for b in range(len(Rr))), sp.Integer(0)) for fr in fun2] for a in al])
        if len(al) == 1 and x != y and not nx.has_path(forest, x, y):
            forest.add_edge(x, y)
            # normalise: leading functional value = 1 -> solve for one symbol with nonzero coefficient
            expr = Mx[0, 0]; b0 = next(b for b in range(len(Rr)) if fun2[0][b] != 0)
            subs[syms[al[0]][b0]] = sp.solve(sp.Eq(expr, 1), syms[al[0]][b0])[0]
        else:
            dets.append(Mx.det())
    eqs = []
    for r in rels:
        p0 = next(iter(r)); X0, Z0 = ix[p0[0][0]], ix[p0[-1][1]]
        tot = [sp.Integer(0)] * len(E.basis[(X0, Z0)])
        for p, c in r.items():
            vec = img[p[0]]; cur = ix[p[0][1]]
            for a in p[1:]:
                nxt = ix[a[1]]; M = E.M[(X0, cur, nxt)]; w = img[a]
                d = len(E.basis[(X0, nxt)])
                vec = [sp.expand(sum((vec[i] * w[j] * M[i][j][c2] for i in range(len(vec)) for j in range(len(w)) if M[i][j][c2]), sp.Integer(0))) for c2 in range(d)]
                cur = nxt
            for c2 in range(len(tot)): tot[c2] += sp.Rational(c.numerator, c.denominator) * vec[c2] if isinstance(c, Fraction) else sp.nsimplify(c) * vec[c2]
        eqs += [sp.expand(e.subs(subs)) for e in tot]
    eqs = [e for e in eqs if e != 0]
    t = sp.Symbol('tRab')
    D = sp.expand(sp.Mul(*[dd.subs(subs) for dd in dets])) if dets else sp.Integer(1)
    if D == 0: return 'NO:det0'
    eqs.append(1 - t * D)
    gens = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
    if not gens: return 'ISO' if all(e != 0 for e in eqs[-1:]) else 'NO:const'
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(tmax)
    try:
        G = sp.groebner(eqs, *gens, order='grevlex', domain='QQ')
        signal.alarm(0)
    except Timeout:
        return witness(eqs, gens, t) if WITNESS else 'UNDECIDED'
    return 'NO:groebner' if list(G.exprs) == [1] else 'ISO'

WITNESS = False
def witness(eqs, gens, t, tries=12, tmax=25):
    """sound ISO certificate only: substitute random integers for some unknowns (not t); a specialised system with a
    common zero over the closure (GB != [1]) gives a zero of the full system.  [1] after specialisation proves nothing."""
    rng = random.Random(7); xs = [g for g in gens if g != t]
    for n in range(tries):
        kk = max(1, int(len(xs) * (0.05 + 0.05 * n)))
        sub = {x: rng.randint(-3, 3) for x in rng.sample(xs, kk)}
        e2 = [sp.expand(e.subs(sub)) for e in eqs]; e2 = [e for e in e2 if e != 0]
        if any(e.is_number for e in e2): continue
        g2 = sorted(set().union(*[e.free_symbols for e in e2]), key=str)
        signal.alarm(tmax)
        try:
            G = sp.groebner(e2, *g2, order='grevlex', domain='QQ'); signal.alarm(0)
        except Timeout:
            continue
        if list(G.exprs) != [1]: return 'ISO(witness)'
    return 'UNDECIDED'

# ---------------- planted wrong children ----------------
def parallel_pair(p, q):
    """p, q differ in exactly one letter and those letters are parallel arrows: then p - l q = u (a - l b) w and the
    substitution a' = a - l b turns it into a monomial, so KILL gives an ISOMORPHIC algebra (not a wrong child)."""
    if len(p) != len(q): return False
    d = [n for n in range(len(p)) if p[n] != q[n]]
    return len(d) == 1 and p[d[0]][:2] == q[d[0]][:2]

def planted(nodes, arrows, rels, how):
    """KILL: replace a binomial p - l q by p (only if the Cartan matrix is unchanged; binomials that are not a
    parallel-letter pair are tried first). SCALE: l -> 2 l.  Returns (rels', tag) or None."""
    order = sorted(range(len(rels)), key=lambda n: len(rels[n]) == 2 and parallel_pair(*sorted(rels[n])))
    for n in order:
        r = rels[n]
        if len(r) != 2: continue
        (p, a), (q, b) = sorted(r.items())
        if how == 'KILL': new = {p: 1}
        else: new = {p: a, q: 2 * b}
        R2 = rels[:n] + [new] + rels[n + 1:]
        if how == 'KILL' and alg_cartan(nodes, arrows, R2) != alg_cartan(nodes, arrows, rels): continue
        return R2, ('par' if parallel_pair(p, q) else 'gen')
    return None

# ---------------- one step ----------------
def step(A, k, wrong=True):
    ch = J.child_of(A, k)
    if ch is None: return None
    nodes, arrows, rels = child_data(ch)
    if not nx.is_directed_acyclic_graph(ch.quiver): return dict(skip='cyclic child')
    Jv = J.J_lib(A, k)
    E = EndT(A, k).table().rad()
    C = E.C; T = E.T
    Hm1 = sum(J.homdim(C, X[1], Yy[1], -1) for X in T for Yy in T)
    Hp1 = sum(J.homdim(C, X[1], Yy[1], 1) for X in T for Yy in T)
    out = dict(J=Jv, Hm1=Hm1, Hp1=Hp1, dimE=sum(E.cartan().values()), narr=len(arrows), nrel=len(rels))
    out['dir'] = iso(E, nodes, arrows, rels)
    out['opp'] = iso(E, *reverse(nodes, arrows, rels))
    if wrong:
        for how in ('KILL', 'SCALE'):
            R2 = planted(nodes, arrows, rels, how)
            out[how] = None if R2 is None else R2[1] + ' ' + iso(E, nodes, arrows, R2[0])
    return out

def summarize(rows):
    tab = collections.Counter()
    for r in rows:
        if r is None or 'skip' in r: tab[('skipped', None if r is None else r['skip'])] += 1; continue
        tab[('J>0' if r['J'] else 'J=0', 'tilt' if r['Hm1'] == 0 and r['Hp1'] == 0 else 'H(%d,%d)' % (r['Hm1'], r['Hp1']),
             'dir ' + r['dir'].split(':')[0], 'opp ' + r['opp'], 'KILL ' + str(r.get('KILL')), 'SCALE ' + str(r.get('SCALE')))] += 1
    for k, v in sorted(tab.items(), key=str): print('  ', v, k)

def controls():
    sq = [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
    SQ = J.build(sq, [{(0, 2): 1, (1, 3): -1}])                     # abd = acd, tilting at 4
    E080 = J.build(sq, [{(0, 2, 4): 1, (1, 3, 4): -1}])              # J != 0 at 4
    A5 = J.build([(1, 2), (2, 3), (3, 4), (4, 5)], [{(0, 1, 2): 1}])
    for name, A, k in (('abd=acd @4', SQ, 4), ('abd=acd @2', SQ, 2), ('E-080 @4 (J!=0)', E080, 4), ('A5/(abc) @2', A5, 2), ('A5/(abc) @3', A5, 3)):
        print(name, mutation.mutationIsPossibleAtVertex(A, k), step(A, k), flush=True)
    # planted: the child of a different vertex
    ch3 = J.child_of(A5, 3); E = EndT(A5, 2).table().rad()
    print('PLANTED child@3 vs End(T@2):', iso(E, *child_data(ch3)))

def lna(nmax):
    from quivermutation import nakayama as nk
    rows = []
    for n in range(4, nmax + 1):
        for A in lna_list(n):
            for B in (A, pathAlgebra.dualPathAlgebra(A)):
                for k in sorted(B.quiver.nodes):
                    if not mutation.mutationIsPossibleAtVertex(B, k): continue
                    rows.append(step(B, k))
    summarize(rows)

def lna_list(n):
    """all LNAs of n vertices with relations given as zero paths (Kupisch series), by brute force over relation sets."""
    out = []
    arrows = [(i, i + 1) for i in range(1, n)]
    # minimal zero relations = intervals [a, b] (arrows a..b-1), length >= 2, no interval containing another
    ivs = [(a, b) for a in range(n - 1) for b in range(a + 2, n)]
    for m in range(0, n):
        for S in itertools.combinations(ivs, m):
            if any(x != y and x[0] <= y[0] and y[1] <= x[1] for x in S for y in S): continue
            out.append(J.build(arrows, [{tuple(range(a, b)): 1} for a, b in S]))
    return out

def rand(N, seed, nmax):
    rng = random.Random(seed); rows = []; t0 = time.time()
    while len(rows) < N:
        al = TH.randomAlg(rng)
        if al.n > nmax: continue
        pa = TH.toLib(al)
        for k in sorted(pa.quiver.nodes):
            if not mutation.mutationIsPossibleAtVertex(pa, k): continue
            t1 = time.time()
            try: rows.append(step(pa, k))
            except AssertionError as e: rows.append(dict(skip='assert ' + str(e)[:30]))
            dt = time.time() - t1
            if dt > 5: print('  slow step %.0f s: n=%d k=%d arrows=%s rels=%s -> %s' % (dt, al.n, k, al.arrows, al.rels, rows[-1]), file=sys.stderr, flush=True)
        if len(rows) % 50 < 7: print('  ... %d steps, %.0f s' % (len(rows), time.time() - t0), file=sys.stderr, flush=True)
    print('time %.0f s' % (time.time() - t0))
    summarize(rows)

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'controls': controls()
    elif mode == 'lna': lna(int(sys.argv[2]))
    elif mode == 'undecided':     # replay the parents whose 'dir' verdict timed out in the rand runs (from the .out files)
        import ast, glob, re
        WITNESS = True
        for f in sorted(glob.glob('workshop/rounds/134/theorist_endt_rand_13%s.out' % (sys.argv[2] if len(sys.argv) > 2 else '?'))):
            for line in open(f):
                if "'dir': 'UNDECIDED'" not in line: continue
                m = re.search(r'n=(\d+) k=(\d+) arrows=(\[.*?\]) rels=(\[.*?\]) -> ', line)
                al = TH.Alg(int(m.group(1)), ast.literal_eval(m.group(3)), ast.literal_eval(m.group(4)))
                pa = TH.toLib(al); k = int(m.group(2)); t1 = time.time()
                r = step(pa, k, wrong=False)
                print(f.split('_')[-1], 'k=%d' % k, r, '%.0f s' % (time.time() - t1), flush=True)
    elif mode == 'rand': rand(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
