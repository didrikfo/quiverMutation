"""Round 053 (toolsmith), T10 item 1: for each of the 13 E-161 edges compute End(T) as a quiver with relations and compare with the next algebra.
Usage (repo root):  .venv/bin/python workshop/rounds/053/toolsmith_endt_run.py [selftest|path13]     (path13 needs /tmp/tsm/c1.pkl from
  DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   (run twice))
"""
import sys, pickle, hashlib
MODE = sys.argv[1] if len(sys.argv) > 1 else 'path13'
sys.argv = ['x', '/tmp/tsm/c1.pkl', '1', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/053/toolsmith_endt.py').read())


def radical_data(TE):
    """rad^m(i,j) as Echelon of coordinate vectors in Hom(i,j) coords; arrows = lifts of rad/rad^2."""
    V = TE.V; dims = {(i, j): TE.dim(i, j) for i in V for j in V}
    for i in V:
        assert dims[(i, i)] == 1, ('End(T_i) not K', i, dims[(i, i)])
    # basis maps for Hom(i,j): TE.basis_maps; coordinate vector of basis t = {t: 1}
    rad = {1: {(i, j): [{t: Fraction(1)} for t in range(dims[(i, j)])] for i in V for j in V if i != j}}
    for m in range(2, len(V) + 1):
        cur = {}
        for i in V:
            for j in V:
                if i == j: continue
                E = Echelon(); vecs = []
                for l in V:
                    if l in (i, j): continue
                    for vec in rad[m - 1].get((i, l), []):
                        f = TE_vec(TE, i, l, vec)
                        for t, g in enumerate(TE.basis_maps(l, j)):
                            h = TE.compose(i, l, j, f, g)
                            cv = TE.coords(i, j, h)
                            if E.add(cv): vecs.append(cv)
                cur[(i, j)] = vecs
        rad[m] = cur
        if not any(cur.values()): break
    return dims, rad


def TE_vec(TE, i, j, vec):
    out = {}
    for t, c in vec.items():
        for kk, x in TE.basis_maps(i, j)[t].items(): out[kk] = out.get(kk, 0) + c * x
    return {k: v for k, v in out.items() if v != 0}


def endquiver(TE):
    dims, rad = radical_data(TE); V = TE.V; arrows = []; layers = {}
    for i in V:
        for j in V:
            if i == j: continue
            E = Echelon()
            for vec in rad.get(2, {}).get((i, j), []): E.add(vec)
            r2 = len(E.rows); picks = []
            for t in range(dims[(i, j)]):
                if E.add({t: Fraction(1)}): picks.append(t)
            for n, t in enumerate(picks): arrows.append((i, j, n, t))
            layers[(i, j)] = [len(rad[m].get((i, j), [])) if m > 1 else dims[(i, j)] for m in sorted(rad)]
    TE.rad = rad
    return dims, arrows, layers


def child_info(c):
    A = Alg(c); V = sorted(c.quiver.nodes)
    d = {(i, j): len(A.basis(i, j)) for i in V for j in V}
    return A, d, ap.arrowsOf(c.quiver), procedure.relationsFrom(c)


def compare(a, v, c, verbose=True):
    """End(T_v(a)) vs c. returns dict of verdicts."""
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); V = TE.V
    A, d, carrows, crels = child_info(c)
    res = {}
    # dims: End Hom(i,j) vs paths i~>j in c (same) or j~>i (opp)
    same = all(dims[(i, j)] == d[(i, j)] for i in V for j in V); opp = all(dims[(i, j)] == d[(j, i)] for i in V for j in V)
    res['dims'] = 'same' if same else ('opp' if opp else 'NO')
    # arrows count per pair
    ec = {}
    for (i, j, n, t) in arrows: ec[(i, j)] = ec.get((i, j), 0) + 1
    cc = {}
    for (s, h, k) in carrows: cc[(s, h)] = cc.get((s, h), 0) + 1
    res['arrows'] = 'same' if ec == cc else ('opp' if ec == {(j, i): n for (i, j), n in cc.items()} else 'NO')
    res['n_arrows'] = (len(arrows), len(carrows))
    mode = res['arrows'] if res['arrows'] != 'NO' else None
    res['rels'] = 'n/a'
    if mode and res['dims'] == res['arrows']:
        res['rels'] = relcheck(TE, arrows, mode, carrows, crels, A)
        if not res['rels'].startswith('all '): res['sym'] = symcheck(TE, arrows, mode, carrows, crels)
    res['layers'] = {k: x for k, x in layers.items() if any(x[1:]) or True}
    res['maxdim'] = max(dims.values())
    return res


def symcheck(TE, arrows, mode, carrows, crels):
    """General decision: do arrow scalars lam_a != 0 and radical-square corrections x_{a,s} exist with all relations of c zero in End(T)?
    Arrow a of c <-> lift b_a = lam_a * (pick) + sum_s x_{a,s} * (rad^2 vector s).  Polynomial system, Groebner basis over Q (with z*prod(lam) = 1).
    Returns 'iso (solution exists over the algebraic closure)' / 'NO label-preserving iso' / 'undecided'."""
    import sympy
    pick = {(i, j): t for (i, j, n, t) in arrows}; cnt = {}
    for (i, j, n, t) in arrows: cnt[(i, j)] = cnt.get((i, j), 0) + 1
    if any(x > 1 for x in cnt.values()): return 'undecided (parallel arrows)'
    syms = []; el = {}
    for ar in carrows:
        s_, h, key = ar; pair = (s_, h) if mode == 'same' else (h, s_)
        lam = sympy.Symbol('l%d_%d' % pair); syms.append(lam)
        vec = {pick[pair]: lam}
        for n_, rv in enumerate(TE.rad.get(2, {}).get(pair, [])):
            x = sympy.Symbol('x%d_%d_%d' % (pair + (n_,))); syms.append(x)
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
    lams = [x for x in syms if str(x).startswith('l')]; z = sympy.Symbol('z')
    eqs.append(z * sympy.prod(lams) - 1)
    if len(eqs) == 1: return 'iso (no equations: all relations vanish identically)'
    G = sympy.groebner(eqs, *(syms + [z]), order='grevlex')
    return 'iso (system with %d unknowns, %d equations consistent)' % (len(syms), len(eqs) - 1) if list(G.exprs) != [1] else 'NO label-preserving iso (1 in the ideal)'


def relcheck(TE, arrows, mode, carrows, crels, A):
    """Send the arrows of c to lifted morphisms and test every relation of c vanishes (scalars: all 1, then diagonal-torus solve if only binomials fail)."""
    # arrow map: c arrow (s,h,key) -> End arrow; multiplicity 1 per pair required
    lift = {}
    for (i, j, n, t) in arrows: lift.setdefault((i, j), []).append(TE.basis_maps(i, j)[t])
    amap = {}
    for (s, h, key) in carrows:
        pair = (s, h) if mode == 'same' else (h, s)
        L = lift[pair]
        if len(L) != 1: return 'parallel-arrows (not decided)'
        amap[(s, h, key)] = (pair, L[0])
    def pathmor(p):
        """c-path (tuple of arrows, first first) -> chain map of End(T), in the category order; returns (i, j, f)."""
        seq = list(p) if mode == 'same' else list(reversed(p))
        (i, j), f = amap[seq[0]][0], amap[seq[0]][1]; cur_i, cur_j = i, j
        for ar in seq[1:]:
            (i2, j2), g = amap[ar]; assert i2 == cur_j; f = TE.compose(cur_i, cur_j, j2, f, g); cur_j = j2
        return cur_i, cur_j, f
    mono = ap.isMonomial(crels) if hasattr(ap, 'isMonomial') else None
    bad = []; nonmono = 0; scal = {}
    for r in crels:
        terms = list(r.items())
        if len(terms) == 1:
            p, cf = terms[0]; i, j, f = pathmor(p)
            if TE.coords(i, j, f) if f else False:
                bad.append(('monomial relation does not vanish', p))
        else:
            nonmono += 1
            vecs = []
            for p, cf in terms:
                i, j, f = pathmor(p); vecs.append((p, cf, i, j, TE.coords(i, j, f) if f else {}))
            scal.setdefault('rel', []).append(vecs)
    if bad: return 'FAIL ' + str(bad[:2])
    if nonmono == 0: return 'all %d relations monomial and vanish in End(T)' % len(crels)
    # non-monomial relations: sum_p cf_p * Lam_p * cv_p = 0 with Lam_p = prod of the arrow scalars on p.  cv_p must be proportional to a common
    # vector (cv_p = rho_p ref); with two nonzero terms this is the binomial Lam_p / Lam_q = -cf_q rho_q / (cf_p rho_p) on the torus (K*)^arrows.
    import sympy
    arr_idx = {a_: n for n, a_ in enumerate(carrows)}; rows = []; targets = []; undecided = []
    for vecs in scal['rel']:
        nz = [(p, cf, cv) for (p, cf, i, j, cv) in vecs if cv]
        if len(nz) < 2: return 'FAIL non-monomial relation with <= 1 nonzero image term %s' % (vecs[0][0],)
        ref = nz[0][2]; k0 = next(iter(ref)); rho = []
        for p, cf, cv in nz:
            if set(cv) != set(ref) or any(cv[k] * ref[k0] != ref[k] * cv[k0] for k in ref): return 'FAIL images not proportional in %s' % (p,)
            rho.append(cv[k0] / ref[k0])
        if len(nz) > 2: undecided.append(len(nz)); continue
        (p, cf, cv), (q, cg, cw) = nz; row = [0] * len(carrows)
        for ar in p: row[arr_idx[ar]] += 1
        for ar in q: row[arr_idx[ar]] -= 1
        rows.append(row); targets.append(-(cg * rho[1]) / (cf * rho[0]))
    if undecided: return 'non-monomial: a relation with %s nonzero terms, torus solve not done' % undecided[:2]
    M = sympy.Matrix(rows); ker = M.T.nullspace()   # y with y M = 0
    for y in ker:
        den = sympy.ilcm(*[sympy.fraction(x)[1] for x in y]); yi = [int(x * den) for x in y]; prod = Fraction(1)
        for e, t in zip(yi, targets): prod *= Fraction(t) ** e
        if prod != 1: return 'FAIL torus inconsistent (product %s)' % prod
    return 'non-monomial relations (%d binomial), scalar torus system consistent (left kernel dim %d): arrow scalars exist' % (len(rows), len(ker))


if MODE == 'selftest':
    # hereditary A3 and a small LNA: End(T_v(a)) vs mutation, all vertices
    cnt = 0
    for lna in nk.LinearNakayamaAlgebra.allOfLength(5)[:12]:
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.quiver.nodes):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                raw = mutation.quiverMutationAtVertex(alg, v)
                c = reduction.reducePathAlgebra(raw)
                if sorted(c.vertices()) != sorted(alg.quiver.nodes): continue
                r = compare(alg, v, c); cnt += 1
                print(v, r['dims'], r['arrows'], r['n_arrows'], r['rels'], 'tiltingPlus', bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v)), 'J', bool(perI(alg, v)))
    print('steps', cnt)
    sys.exit()

def edges13():
    recs = [r for r in pickle.load(open('/tmp/tsm/c1.pkl', 'rb')) if r['kind'] == 'fail']
    out = []; x0 = recs[14]['childObj']
    for name, x, mv in [('child', x0, 'F1 F3 F1 F5 F4 R7 R1'.split()), ('lna9', classes[base][9], 'F7 R2 R1 F2 F7 R3'.split())]:
        for s_ in mv:
            kind, v = s_[0], int(s_[1:]); a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
            c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v)); out.append((name, s_, a, v, c))
            x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return out

if MODE == 'perturb':
    # power control of the iso decision: (A) one binomial relation of c with a coefficient doubled, (B) one binomial replaced by its first term
    # (a monomial the End(T) does not satisfy).  Expected: B -> NO always; A -> NO when the torus system has a cycle through it.
    tot = dict(A_no=0, A_iso=0, B_no=0, B_iso=0)
    for name, s_, a, v, c in edges13():
        TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); carrows = ap.arrowsOf(c.quiver); crels = procedure.relationsFrom(c)
        res = []
        for n_, r in enumerate(crels):
            if len(r) < 2: continue
            pa = list(crels); k0 = next(iter(r)); rr = dict(r); rr[k0] = rr[k0] * 2; pa[n_] = rr
            pb = list(crels); pb[n_] = {k0: r[k0]}
            va = symcheck(TE, arrows, 'same', carrows, pa); vb = symcheck(TE, arrows, 'same', carrows, pb)
            tot['A_no' if va.startswith('NO') else 'A_iso'] += 1; tot['B_no' if vb.startswith('NO') else 'B_iso'] += 1
            res.append((n_, va[:2], vb[:2]))
        print(name, s_, res, flush=True)
    print('PERTURB', tot); sys.exit()

if MODE == 'tilt':   # Hom(T,T[m]) for m = -1, 0, 1 with the skeptic's code, tie to End dims
    _st = open('workshop/rounds/050/skeptic_tilt.py').read(); P = 2147483647; exec(_st[_st.index('def rank'):])
    for name, s_, a, v, c in edges13():
        Vs, res = stepTest(a, v, -1); TE = TiltEnd(a, v)
        H = [[TE.dim(i, j) for j in Vs] for i in Vs]
        print(name, s_, 'Hom(T,T[-1])', sum(map(sum, res[-1])), 'Hom(T,T[1])', sum(map(sum, res[1])), 'H == own End dims', H == res[0], flush=True)
    sys.exit()

def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]
recs = [r for r in pickle.load(open('/tmp/tsm/c1.pkl', 'rb')) if r['kind'] == 'fail']
if MODE == 'fail':     # negative control: the failing J != 0 key-keeping steps (parent -> child), End(T_v(parent)) vs child
    for n_, r in enumerate(recs):
        a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
        rr = compare(a, v, c); print('fail', n_, 'v', v, 'dims', rr['dims'], 'arrows', rr['arrows'], rr['n_arrows'], '|', rr['rels'], flush=True)
    sys.exit()
CS = 'F1 F3 F1 F5 F4 R7 R1'.split(); TS = 'F7 R2 R1 F2 F7 R3'.split()
x0 = recs[14]['childObj']; print('child key', tag(x0))
for name, x, mv in [('child', x0, CS), ('lna9', classes[base][9], TS)]:
    for s_ in mv:
        kind, v = s_[0], int(s_[1:])
        a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        raw = mutation.quiverMutationAtVertex(a, v); c = reduction.reducePathAlgebra(raw)
        r = compare(a, v, c)
        print(name, s_, 'dims', r['dims'], 'arrows', r['arrows'], r['n_arrows'], 'max dim Hom', r['maxdim'], '|', r['rels'][:60], '|', r.get('sym', ''), flush=True)
        x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
