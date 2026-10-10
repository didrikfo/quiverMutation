"""Round 057 (toolsmith), T10 item 1: End(T) as a quiver with relations, now with MATRIX-VALUED arrow identification (parallel arrows).
Library for toolsmith_endt2_run.py.  Reuses the round-053 code (workshop/rounds/053/toolsmith_endt.py, ..._run.py: TiltEnd, endquiver, child_info).
symcheck2: for every pair (i,j) with m arrows of c (parallel if m > 1) the arrows b_1..b_m of c are sent to
   b_k -> sum_l M[k][l] * pick_l + sum_s x[k][s] * (rad^2 basis vector s)         (pick_l: lifts of a basis of rad/rad^2 of End(T)(i,j))
with M an m x m matrix of unknowns, det M != 0 enforced by z * prod_pairs det(M) = 1; every relation of c must vanish in End(T).
A Groebner basis != [1] gives a solution over the algebraic closure => K Q_c/I_c -> End(T) is onto (M invertible, rad^2 corrections) and between
equal finite dimensions (checked: dims) => isomorphism, label-preserving on vertices; if the basis is [1] there is NO iso that preserves vertex labels.
For m = 1 on every pair this is the round-053 symcheck."""
import sys, os
CLS = int(os.environ.get('CLS', '1'))
PK = os.environ.get('PK', '/tmp/tsm/c%d.pkl' % CLS)
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/053/toolsmith_endt.py').read())
_r = open('workshop/rounds/053/toolsmith_endt_run.py').read()
exec(_r[_r.index('def radical_data'):_r.index('def compare(')])   # radical_data, TE_vec, endquiver, child_info
import sympy


def symcheck2(TE, arrows, mode, carrows, crels, timeit=None):
    pickl = {}
    for (i, j, n, t) in arrows: pickl.setdefault((i, j), []).append(t)
    cgroups = {}
    for ar in carrows:
        s_, h, key = ar; pair = (s_, h) if mode == 'same' else (h, s_); cgroups.setdefault(pair, []).append(ar)
    if {p: len(v) for p, v in pickl.items()} != {p: len(v) for p, v in cgroups.items()}: return 'NO (arrow counts differ)'
    syms = []; el = {}; dets = []
    for pair, ars in sorted(cgroups.items()):
        m = len(ars); M = sympy.Matrix(m, m, lambda k, l: sympy.Symbol('m%d_%d_%d_%d' % (pair + (k, l)))); syms += list(M)
        dets.append(M.det())
        for k, ar in enumerate(ars):
            vec = {}
            for l in range(m): vec[pickl[pair][l]] = M[k, l]
            for n_, rv in enumerate(TE.rad.get(2, {}).get(pair, [])):
                x = sympy.Symbol('x%d_%d_%d_%d' % (pair + (k, n_))); syms.append(x)
                for t, c in rv.items(): vec[t] = vec.get(t, 0) + sympy.Rational(c.numerator, c.denominator) * x
            el[ar] = (pair, vec)
    cache = {}
    def mult(i, l, j, t1, t2):
        if (i, l, j, t1, t2) not in cache:
            h = TE.compose(i, l, j, TE.basis_maps(i, l)[t1], TE.basis_maps(l, j)[t2])
            cache[(i, l, j, t1, t2)] = TE.coords(i, j, h) if h else {}
        return cache[(i, l, j, t1, t2)]
    def pathval(p):
        seq = list(p) if mode == 'same' else list(reversed(p))
        (i, l), f = el[seq[0]]
        for ar in seq[1:]:
            (l2, j), g = el[ar]; assert l2 == l; out = {}
            for t1, a1 in f.items():
                for t2, a2 in g.items():
                    for t, c in mult(i, l, j, t1, t2).items(): out[t] = out.get(t, 0) + sympy.Rational(c.numerator, c.denominator) * a1 * a2
            f = {t: sympy.expand(v) for t, v in out.items()}; l = j
        return (i, l), f
    eqs = []
    for r in crels:
        tot = {}
        for p, cf in r.items():
            _, f = pathval(p)
            for t, v in f.items(): tot[t] = tot.get(t, 0) + sympy.Rational(cf.numerator, cf.denominator) * v
        eqs += [sympy.expand(v) for v in tot.values() if sympy.expand(v) != 0]
    nrel = len(eqs); zs = [sympy.Symbol('z%d' % n_) for n_ in range(len(dets))]; eqs += [sympy.expand(zz * dd - 1) for zz, dd in zip(zs, dets)]   # one Rabinowitsch variable per block (cheaper than z * product)
    nm = sum(1 for g in cgroups.values() if len(g) > 1)
    if nrel == 0: return 'iso (no equations; %d parallel pairs)' % nm
    G = sympy.groebner(eqs, *(syms + zs), order='grevlex')
    if list(G.exprs) == [1]: return 'NO label-preserving iso (1 in the ideal; %d parallel pairs)' % nm
    return 'iso (%d unknowns, %d equations consistent; %d parallel pairs)' % (len(syms), len(eqs) - 1, nm)


def compare2(a, v, c):
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); V = TE.V
    A, d, carrows, crels = child_info(c); res = {}
    same = all(dims[(i, j)] == d[(i, j)] for i in V for j in V); opp = all(dims[(i, j)] == d[(j, i)] for i in V for j in V)
    res['dims'] = 'same' if same else ('opp' if opp else 'NO')
    ec = {}
    for (i, j, n, t) in arrows: ec[(i, j)] = ec.get((i, j), 0) + 1
    cc = {}
    for (s, h, k) in carrows: cc[(s, h)] = cc.get((s, h), 0) + 1
    res['arrows'] = 'same' if ec == cc else ('opp' if ec == {(j, i): n for (i, j), n in cc.items()} else 'NO')
    res['n_arrows'] = (len(arrows), len(carrows)); res['par'] = max(cc.values()) if cc else 0; res['maxdim'] = max(dims.values())
    res['verdict'] = 'NO (dims/arrows)' if 'NO' in (res['dims'], res['arrows']) or res['dims'] != res['arrows'] else symcheck2(TE, arrows, res['arrows'], carrows, crels)
    return res


def edges_of(start, moves_, tags=''):
    """moves_: ['F1', 'R7', ...] applied forward from `start` (an algebra); yields (label, a, v, c) for each edge (a = algebra or its dual, v vertex, c the reduced mutation)."""
    x = start
    for s_ in moves_:
        kind, v = s_[0], int(s_[1:]); a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v)); yield (tags + s_, a, v, c)
        x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
