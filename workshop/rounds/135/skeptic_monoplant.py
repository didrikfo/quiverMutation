"""Round 135, skeptic, T5 / H-023: planted controls for E-345's child = End(T) isomorphism test (`iso` of
rounds/134/theorist_endt.py, used unchanged).

Two families of WRONG children with the child's quiver and the child's (= End(T)'s) labelled Cartan matrix:

 mono    : at every gate-admitted step with an acyclic, MONOMIAL, library child, every other monomial ideal on the
           child's quiver (antichain of zero paths of length >= 2) with the same labelled Cartan matrix.
           Truth = NO when the quiver has no parallel arrows (proof: a label-preserving iso phi fixes e_x up to
           conjugation and sends an arrow a to lambda_a a + (longer paths), lambda_a != 0 because the arrow slot is
           1-dimensional mod rad^2; so a path p of length m goes to lambda p + (paths of length > m); in a monomial
           algebra the nonzero paths are a basis, hence p is zero in one algebra iff it is zero in the other).
           With parallel arrows: the label swap a0 <-> a1 of the child's relations is a TRUE-ISO control; other
           monomial ideals there have no truth (Groebner verdict only, reported apart).
 modulus : E-186-type.  Parent M(p) on 7 vertices: w -> x, a0, a1: x -> y, c_j: y -> z_j (j = 1..4), relations
           (s_j a0 + t_j a1) c_j = 0 for four points p_j = [s_j : t_j] of P^1.  Mutation at the source w.  The child's
           relations among {x, y, z_j} are the same four point relations (checked at run time); a label-preserving iso of
           children restricts to the corner algebra e C e (e = all vertices but w; checked: w is a sink or source of the
           child quiver, so no path through w returns), where it acts on span{a0, a1} by GL2 and on each c_j by a scalar.
           So truth: ISO iff the ordered 4-tuples of points are PGL2-equivalent (distinct points: same cross-ratio;
           coincidences are invariant).  Planted children: the child with its point relations replaced by another tuple.
usage (repo root):  .venv/bin/python -W ignore workshop/rounds/135/skeptic_monoplant.py modulus
                    .venv/bin/python -W ignore workshop/rounds/135/skeptic_monoplant.py mono lna 6
                    .venv/bin/python -W ignore workshop/rounds/135/skeptic_monoplant.py mono rand N SEED
"""
import sys, collections, itertools, random, time
from fractions import Fraction
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/134')
import warnings; warnings.simplefilter('ignore')
import networkx as nx
import theorist_endt as T        # iso, EndT, child_data, alg_cartan, lna_list, J (jsilt), TH (thin)
from quivermutation import mutation, pathAlgebra
J, TH = T.J, T.TH

# ---------------- monomial planted children ----------------
def allpaths(arrows, nodes, minlen=2):
    out = []
    for x in nodes:
        for y in nodes:
            if x != y: out += [p for p in TH.paths(arrows, x, y) if len(p) >= minlen]
    return out

def contains(p, q):
    """q is a subpath (consecutive letters) of p"""
    return any(p[i:i + len(q)] == q for i in range(len(p) - len(q) + 1))

def mono_cartan(nodes, arrows, S, allp):
    c = collections.Counter()
    for p in allp:
        if not any(contains(p, q) for q in S): c[(p[0][0], p[-1][1])] += 1
    for x in nodes:
        c[(x, x)] = 1
        for y in nodes:
            if x != y: c[(x, y)] += sum(1 for a in arrows if a[0] == x and a[1] == y) * 0   # arrows are length-1 paths
    for a in arrows: c[(a[0], a[1])] += 1
    return {(x, y): c.get((x, y), 0) for x in nodes for y in nodes}

def mono_alternatives(nodes, arrows, rels, cap=40):
    """all antichains S of paths (len >= 2), S != child's zero paths, with the child's Cartan matrix (up to cap)."""
    allp = allpaths(arrows, nodes)
    own = frozenset(next(iter(r)) for r in rels)
    target = mono_cartan(nodes, arrows, own, allp)
    # deficit per slot: number of zero paths per slot must equal the child's
    need = {xy: sum(1 for p in allp if (p[0][0], p[-1][1]) == xy) + sum(1 for a in arrows if (a[0], a[1]) == xy) - target[xy]
            for xy in target if xy[0] != xy[1]}
    out = []
    cand = sorted(allp, key=len)
    def rec(i, S):
        if len(out) >= cap: return
        if i == len(cand):
            if frozenset(S) != own and mono_cartan(nodes, arrows, S, allp) == target: out.append(list(S))
            return
        p = cand[i]
        rec(i + 1, S)
        if not any(contains(p, q) for q in S):
            # prune: adding p must not over-kill a slot (killing p kills all its superpaths)
            S2 = S + [p]
            c = mono_cartan(nodes, arrows, S2, allp)
            if all(c[xy] >= target[xy] for xy in target): rec(i + 1, S2)
    rec(0, [])
    return out

EX = []
def mono_truth(arrows, I, I2):
    """label-preserving iso of two admissible monomial algebras kQ/I, kQ/I2 on one acyclic quiver: ISO iff a
    permutation of parallel arrows inside their slots carries the minimal zero paths of I onto those of I2.
    Proof: an iso induces a graded label-preserving iso of gr = the algebras themselves (monomial = path-length graded),
    i.e. g in G = prod GL(arrow slot) with g I = I2.  I, I2 are fixed by the diagonal torus T; if g I = I2 then
    g^-1 T g and T are maximal tori of Stab(I)^0, so g^-1 T g = h T h^-1 with h in Stab(I)^0, gh normalises T, and
    (gh) I = I2 with gh = permutation x torus; the torus fixes monomial ideals."""
    slots = collections.defaultdict(list)
    for a in arrows: slots[(a[0], a[1])].append(a)
    groups = [al for al in slots.values() if len(al) > 1]
    for perms in itertools.product(*[list(itertools.permutations(al)) for al in groups]):
        m = {}
        for al, pm in zip(groups, perms): m.update(dict(zip(al, pm)))
        if frozenset(tuple(m.get(a, a) for a in p) for p in I) == I2: return 'ISO'
    return 'NO'

def has_parallel(arrows):
    return any(v > 1 for v in collections.Counter((a[0], a[1]) for a in arrows).values())

def swap_controls(arrows, rels):
    """relabel one parallel pair a0 <-> a1 in the child's relations (an automorphism of the quiver): truth ISO."""
    out = []
    slots = collections.defaultdict(list)
    for a in arrows: slots[(a[0], a[1])].append(a)
    for al in slots.values():
        if len(al) < 2: continue
        m = {al[0]: al[1], al[1]: al[0]}
        R2 = [{tuple(m.get(a, a) for a in p): c for p, c in r.items()} for r in rels]
        if sorted(map(sorted, (r.items() for r in R2))) != sorted(map(sorted, (r.items() for r in rels))): out.append(R2)
    return out

def steps_lna(nmax):
    for n in range(4, nmax + 1):
        for A in T.lna_list(n):
            for B in (A, pathAlgebra.dualPathAlgebra(A)):
                for k in sorted(B.quiver.nodes):
                    if mutation.mutationIsPossibleAtVertex(B, k): yield ('lna%d' % n, B, k)

def steps_rand(N, seed, nmax=6):
    rng = random.Random(seed); m = 0
    while m < N:
        al = TH.randomAlg(rng)
        if al.n > nmax: continue
        pa = TH.toLib(al)
        for k in sorted(pa.quiver.nodes):
            if mutation.mutationIsPossibleAtVertex(pa, k): m += 1; yield ('rand', pa, k)

def mono(gen):
    tab = collections.Counter(); t0 = time.time(); bad = []
    for name, A, k in gen:
        ch = J.child_of(A, k)
        if ch is None: tab['no child'] += 1; continue
        if not nx.is_directed_acyclic_graph(ch.quiver): tab['cyclic child'] += 1; continue
        nodes, arrows, rels = T.child_data(ch)
        if any(len(r) != 1 for r in rels): tab['child not monomial'] += 1; continue
        tab['monomial child'] += 1
        alts = mono_alternatives(nodes, arrows, rels)
        swaps = swap_controls(arrows, rels) if has_parallel(arrows) else []
        if not alts and not swaps: tab['  no same-Cartan alternative'] += 1; continue
        E = T.EndT(A, k).table().rad()
        own = T.iso(E, nodes, arrows, rels)
        tab['  child itself: ' + own] += 1
        par = has_parallel(arrows)
        ownset = frozenset(next(iter(r)) for r in rels)
        for S in alts:
            v = T.iso(E, nodes, arrows, [{p: 1} for p in S])
            truth = mono_truth(arrows, ownset, frozenset(S))
            key = 'planted mono, %s arrows, truth %s: %s' % ('parallel' if par else 'no parallel', truth, v)
            tab[key] += 1
            if v.startswith('NO'): tab['  stage of NO: ' + v] += 1
            if v.split(':')[0] != truth: bad.append((name, k, arrows, rels, S, v, truth))
            if truth == 'NO' and v.startswith('NO') and len(EX) < 3 and (par or not EX): EX.append((name, k, arrows, rels, S, v))
        for R2 in swaps:
            v = T.iso(E, nodes, arrows, R2)
            tab['swap control, truth ISO: ' + v] += 1
            if v != 'ISO': bad.append((name, k, arrows, rels, R2, v))
    print('time %.0f s' % (time.time() - t0))
    for key, v in sorted(tab.items()): print('  %5d  %s' % (v, key))
    for b in bad: print('DISAGREE', b)
    for e in EX: print('EXAMPLE truth-NO planted:', e)

# ---------------- E-186-type modulus child ----------------
W, X, Y, Z = 0, 1, 2, (3, 4, 5, 6)
def modulus_parent(pts):
    arrows = [(W, X, 0), (X, Y, 0), (X, Y, 1)] + [(Y, z, 0) for z in Z]
    a0, a1 = (X, Y, 0), (X, Y, 1)
    rels = []
    for (s, t), z in zip(pts, Z):
        c = (Y, z, 0); r = {}
        if s: r[(a0, c)] = Fraction(s)
        if t: r[(a1, c)] = Fraction(t)
        rels.append(r)
    return TH.toLib(TH.Alg(7, arrows, rels))          # library labels = 1 + these

def cross(pts):
    """truth class of an ordered 4-tuple of points of P^1: coincidence pattern, and the cross-ratio if distinct."""
    def eq(p, q): return p[0] * q[1] == p[1] * q[0]
    pat = tuple(min(j for j in range(4) if eq(pts[i], pts[j])) for i in range(4))
    if pat != (0, 1, 2, 3): return ('pattern', pat)
    def det(p, q): return Fraction(p[0] * q[1] - p[1] * q[0])
    p1, p2, p3, p4 = pts
    return ('cr', det(p1, p3) * det(p2, p4) / (det(p1, p4) * det(p2, p3)))

def point_rels(rels, a0, a1, cs, pts):
    """replace each relation supported on {a0 c_j, a1 c_j} by the point relation of pts[j]"""
    out = []; done = set()
    for r in rels:
        supp = set(r)
        j = next((j for j, c in enumerate(cs) if supp <= {(a0, c), (a1, c)}), None)
        if j is None: out.append(r); continue
        s, t = pts[j]; nr = {}
        if s: nr[(a0, cs[j])] = Fraction(s)
        if t: nr[(a1, cs[j])] = Fraction(t)
        out.append(nr); done.add(j)
    assert done == set(range(4)), 'child lost a point relation'
    return out

def child_points(rels, a0, a1, cs):
    pts = [None] * 4
    for r in rels:
        for j, c in enumerate(cs):
            if set(r) <= {(a0, c), (a1, c)}:
                s, t = r.get((a0, c), 0), r.get((a1, c), 0); pts[j] = (s, t)
    return pts

def modulus():
    base = [(1, 0), (0, 1), (1, 1), (1, 2)]          # points [s:t] meaning (s a0 + t a1) c_j = 0
    variants = {
        'same tuple (identity)': base,
        'Moebius image, same cross-ratio': None,       # filled below
        'lambda 2 -> 3': [(1, 0), (0, 1), (1, 1), (1, 3)],
        'lambda 2 -> -1': [(1, 0), (0, 1), (1, 1), (1, -1)],
        'lambda 2 -> 1/2': [(1, 0), (0, 1), (1, 1), (2, 1)],
        'p3, p4 swapped': [(1, 0), (0, 1), (1, 2), (1, 1)],
        'p4 = p3 (coincidence)': [(1, 0), (0, 1), (1, 1), (1, 1)],
        'MONOMIAL 0,inf,0,inf': [(1, 0), (0, 1), (1, 0), (0, 1)],
        'MONOMIAL 0,inf,inf,0': [(1, 0), (0, 1), (0, 1), (1, 0)],
        'MONOMIAL 0,0,0,0': [(1, 0), (1, 0), (1, 0), (1, 0)],
    }
    # Moebius map u -> (2u + 1) / (u + 3) on [s:t] with u = t/s... act on (s, t) by a matrix in GL2
    g = ((2, 1), (1, 3))
    variants['Moebius image, same cross-ratio'] = [(g[0][0] * s + g[0][1] * t, g[1][0] * s + g[1][1] * t) for s, t in base]
    pa = modulus_parent(base)
    for k in sorted(pa.quiver.nodes):
        print('vertex', k, 'gate', mutation.mutationIsPossibleAtVertex(pa, k))
    k = W + 1
    ch = J.child_of(pa, k)
    nodes, arrows, rels = T.child_data(ch)
    print('child arrows', arrows); print('child rels', rels)
    # identify a0, a1 (the parallel pair x -> y) and c_j (y -> z_j) in the child, and check the corner claim
    x, y = X + 1, Y + 1
    par = [a for a in arrows if (a[0], a[1]) == (x, y)]
    cs = [next(a for a in arrows if (a[0], a[1]) == (y, z + 1)) for z in Z]
    assert len(par) == 2
    w = W + 1
    ins = [a for a in arrows if a[1] == w]; outs = [a for a in arrows if a[0] == w]
    print('w in-arrows', ins, 'out-arrows', outs, '-> w is a', 'sink' if not outs else 'source' if not ins else 'MIDDLE (corner argument fails)')
    assert not (ins and outs)
    other = [r for r in rels if not any(set(r) <= {(par[0], c), (par[1], c)} for c in cs)]
    print('relations besides the point relations:', other, '(must avoid the corner: each involves w)')
    assert all(any(w in (a[0], a[1]) for p in r for a in p) for r in other)
    cpts = child_points(rels, par[0], par[1], cs)
    print('child points', cpts, 'class', cross(cpts))
    E = T.EndT(pa, k).table().rad()
    print('J', J.J_lib(pa, k), 'child itself vs End(T):', T.iso(E, nodes, arrows, rels))
    # the planted children are written in the child's own coordinates: the child's point tuple is mapped as base is
    # mapped to the variant by position (variants are given relative to base; child points are base up to GL2 and
    # rescaling, checked by the cross-ratio line above)
    assert cross(cpts) == cross(base)
    # map variant tuples into the child's coordinates: find h in GL2 with h(base_j) ~ cpts_j for j = 1, 2, 3
    h = gl2_fit(base[:3], cpts[:3])
    agree = 0; tot = 0
    for name, v in variants.items():
        vv = [apply(h, p) for p in v]
        truth = 'ISO' if cross(vv) == cross(cpts) else 'NO'
        R2 = [{p: c for p, c in r.items()} for r in point_rels(rels, par[0], par[1], cs, vv)]
        t1 = time.time(); verdict = T.iso(E, nodes, arrows, R2)
        ok = verdict.split(':')[0] == truth; agree += ok; tot += 1
        print('  %-34s points %-28s truth %-3s test %-12s %s  %.1f s' % (name, [(int(a) if a.denominator == 1 else str(a), int(b) if b.denominator == 1 else str(b)) for a, b in vv], truth, verdict, 'ok' if ok else 'DISAGREE', time.time() - t1))
    print('agree %d of %d (named variants)' % (agree, tot))
    # random batch: 40 random GL2 images of the child's tuple (truth ISO) and 40 random integer tuples (truth by class)
    rng = random.Random(135); tab = collections.Counter()
    for n in range(80):
        if n < 40:
            while True:
                g = ((rng.randint(-4, 4), rng.randint(-4, 4)), (rng.randint(-4, 4), rng.randint(-4, 4)))
                if g[0][0] * g[1][1] - g[0][1] * g[1][0]: break
            vv = [apply(g, p) for p in cpts]
        else:
            vv = []
            while len(vv) < 4:
                p = (rng.randint(-3, 3), rng.randint(-3, 3))
                if p != (0, 0): vv.append(p)
        truth = 'ISO' if cross(vv) == cross(cpts) else 'NO'
        verdict = T.iso(E, nodes, arrows, point_rels(rels, par[0], par[1], cs, vv))
        tab[(truth, verdict, 'monomial' if all(0 in p for p in vv) else '')] += 1
    for key, v in sorted(tab.items()): print('  random batch: truth %s, test %s %s: %d' % (key + (v,)))

def apply(h, p):
    return (h[0][0] * p[0] + h[0][1] * p[1], h[1][0] * p[0] + h[1][1] * p[1])

def gl2_fit(src, dst):
    """h in GL2(Q) with h(src_j) proportional to dst_j, j = 0, 1, 2 (src, dst three distinct points)."""
    import sympy as sp
    a, b, c, d, l0, l1, l2 = sp.symbols('a b c d l0 l1 l2')
    eqs = []
    for (s, t), (u, v), l in zip(src, dst, (l0, l1, l2)):
        eqs += [a * s + b * t - l * u, c * s + d * t - l * v]
    sol = sp.solve(eqs + [l0 - 1], [a, b, c, d, l0, l1, l2], dict=True)[0]
    return ((Fraction(str(sol[a])), Fraction(str(sol[b]))), (Fraction(str(sol[c])), Fraction(str(sol[d]))))

if __name__ == '__main__':
    MYARGV = list(sys.argv)
    mode = MYARGV[1]
    if mode == 'modulus': modulus()
    elif mode == 'mono':
        if MYARGV[2] == 'lna': mono(steps_lna(int(MYARGV[3])))
        else: mono(steps_rand(int(MYARGV[3]), int(MYARGV[4])))
