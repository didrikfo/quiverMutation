"""Round 123 (theorist), T5 / H-023 / H-031: the round-122 hypotheses (LF), (H1), (H3) on REAL Construction-2.2 parents.
Parent = a minimal model of a small LNA (line 1 -> ... -> n, zero relations), built here cell by cell (Tate-style, exact
rationals): letters of degree m+1 in cell (s, t) kill H_m of the words of shorter letters; the result is checked to have
H = A (dim 1 in degree 0 where the path is nonzero in A, 0 elsewhere).  Weight = path length (arrow weight 1).
Def. 3.2 kM at vertex i (A = the degree-0 letters out of i) from workshop/rounds/069/scholar_oppermann_sign.py, signs
last = -1 (E-202), through workshop/rounds/106/theorist_loopfilt.setup; checks reused from workshop/rounds/122/theorist_h031.py.
Both A and A^op are built at every vertex, so whichever side Def. 3.2 is of the repo's F / R (E-252), the step is covered.
Parts: [MM] minimal model H = A;  [D2] d^2 = 0 on kM;  [LF] positive cycles (Floyd-Warshall);  [H1] linear part of d;
[PI] pi: kM -> T(H) dg map, d_H^2 = 0;  [SH] Lemma 3 shape (g found by recursion over d_H, verified by h031.shape);
[ENG] H(I) cell by cell mod 1000003 (copy of lf.run's cell loop, vertices 1..n), I = (Theta, d Theta).
(H1) is printed strict (E-319) and weak (H1w: only Theta letters need pivots; what Lemmas 1-2 use).  [LF] also on a
two-A-arrow parent with unequal weights (must say NO).  [TREE] a non-LNA gate-admitted J = 0 step where strict (H1) fails.
Controls: [MM] a planted wrong d (drop a term) must fail H = A or d^2 = 0; [SH] planted g; [ENG] the closed ideal (a*)
in place of I must give H != 0.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/123/theorist_lnaparent.py [--rows n:row,...]"""
import sys, importlib.util, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '.')
def load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
h = load('h031', 'workshop/rounds/122/theorist_h031.py'); lf = h.lf; opp = lf.opp

# ---------------- exact linear algebra
def rref_add(basis, v):
    """basis: dict pivot -> row (dict); reduce v; if independent add it and return True"""
    v = {k: Fr(c) for k, c in v.items() if c}
    for p, row in basis.items():
        if p in v:
            f = v[p]
            for k, c in row.items():
                v[k] = v.get(k, 0) - f * c
                if v[k] == 0: del v[k]
    if not v: return False
    p = min(v, key=repr); f = v[p]; v = {k: c / f for k, c in v.items()}
    for q, row in basis.items():
        if p in row:
            g = row[p]
            for k, c in v.items():
                row[k] = row.get(k, 0) - g * c
                if row[k] == 0: del row[k]
    basis[p] = v; return True

def kernel(cols, vecs):
    """vecs: list of dicts (images of basis elements cols); return basis of kernel as dicts col -> coeff"""
    rows = {}; comb = {}   # echelon of images with their combinations
    ker = []
    for c, v in zip(cols, vecs):
        v = {k: Fr(x) for k, x in v.items() if x}; cv = {c: Fr(1)}
        while v:
            p = min(v, key=repr)
            if p not in rows: break
            f = v[p]
            for k, x in rows[p].items():
                v[k] = v.get(k, 0) - f * x
                if v[k] == 0: del v[k]
            for k, x in comb[p].items():
                cv[k] = cv.get(k, 0) - f * x
                if cv[k] == 0: del cv[k]
        if not v: ker.append(cv)
        else:
            p = min(v, key=repr); f = v[p]
            rows[p] = {k: x / f for k, x in v.items()}; comb[p] = {k: x / f for k, x in cv.items()}
    return ker

def integral(v):
    from math import lcm, gcd
    L = 1
    for c in v.values(): L = lcm(L, Fr(c).denominator)
    w = {k: int(c * L) for k, c in v.items()}
    g = 0
    for c in w.values(): g = gcd(g, c)
    return {k: c // g for k, c in w.items()}

# ---------------- minimal model of an LNA
def lna(n, row, op=False):
    rels = [(s, s + l) for s, l in enumerate(row, start=1) if l]       # zero relation = path s -> s+l
    return n, rels, op

def minimal_model(n, rels, op, plant=False):
    arr = {}; d = {}
    E = lambda s, t: (t, s) if op else (s, t)                             # op: every arrow reversed
    def nonzero(s, t):   # path between line vertices s <= t nonzero in A
        return not any(s <= a and b <= t for a, b in rels)
    q = opp.Q(arr, [], 0); q.d = d
    def words(s, t):
        out = []
        def rec(v, w):
            if v == t and len(w) >= 2: out.append(tuple(w))
            for L, (a, b, _) in arr.items():
                if a == v and (b - a) * (t - s) > 0 and abs(b - s) <= abs(t - s): rec(b, w + [L])
        rec(s, []); return out
    cnt = 0
    for wdt in range(1, n):
        for lo in range(1, n - wdt + 1):
            s, t = E(lo, lo + wdt)
            if wdt == 1:
                arr['x%d' % lo] = (s, t, 0); continue
            ws = words(s, t); bydeg = {}
            for w in ws: bydeg.setdefault(q.deg(w), []).append(w)
            for m in sorted(bydeg):
                Z = kernel(bydeg[m], [q.dpath(w) for w in bydeg[m]]) if m > 0 else [{w: 1} for w in bydeg[m]]
                B = {}
                for w in bydeg.get(m + 1, []): rref_add(B, q.dpath(w))
                reps = []
                for z in Z:
                    if rref_add(B, z): reps.append(z)
                if m == 0 and nonzero(lo, lo + wdt): reps = reps[1:]          # keep the one class of the nonzero path
                for z in reps:
                    cnt += 1; nm_ = 'r%d' % cnt if m == 0 else ('c%d_%d' % (m + 1, cnt))
                    arr[nm_] = (s, t, m + 1); d[nm_] = integral(z)
    if plant:   # control: drop one term of the first letter with >= 2 terms in d
        L = next(L for L in d if len(d[L]) >= 2); d[L].pop(next(iter(d[L])))
    return arr, d

def check_mm(n, rels, op, arr, d):
    """H(parent) = A in every cell and d^2 = 0 on letters"""
    q = opp.Q(arr, [], 0); q.d = d
    bad = [L for L in d if q.dlin(d[L])]
    E = lambda s, t: (t, s) if op else (s, t)
    for wdt in range(1, n):
        for lo in range(1, n - wdt + 1):
            s, t = E(lo, lo + wdt)
            ws = []
            def rec(v, w):
                if v == t and w: ws.append(tuple(w))
                for L, (a, b, _) in arr.items():
                    if a == v and (b - a) * (t - s) > 0 and abs(b - s) <= abs(t - s): rec(b, w + [L])
            rec(s, [])
            bydeg = {}
            for w in ws: bydeg.setdefault(q.deg(w), []).append(w)
            for m in sorted(bydeg):
                z = len(kernel(bydeg[m], [q.dpath(w) for w in bydeg[m]])) if m > 0 else len(bydeg[m])
                B = {}
                for w in bydeg.get(m + 1, []): rref_add(B, q.dpath(w))
                hm = z - len(B)
                want = 1 if (m == 0 and not any(lo <= a and b <= lo + wdt for a, b in rels)) else 0
                if hm != want: bad.append(((s, t), m, hm, want))
    return bad

# ---------------- kM at vertex i
def kM(arr, d, i):
    A = sorted(L for L, (s, t, g) in arr.items() if s == i and g == 0)
    q = lf.mk(dict(arr), A, {k: dict(v) for k, v in d.items()}, {L: abs(t - s) for L, (s, t, g) in arr.items()})
    q.i = i
    m, letters, wt = lf.setup(q)
    dval = {L: {w: k for w, k in m.partial({(L,): 1}).items() if k} for L in letters}
    theta = {L for L in letters if L[0] == 'm' and L[3] is not None}
    return q, m, letters, wt, dval, theta

def g_recursive(m, H, dH, theta, base_star, base0):
    """g(x) = 1 + max over words of d_H x of sum g; theta: g = top; letters with d_H = 0: base values"""
    G = {}; stack = set()
    def val(L):
        if L in G: return G[L]
        if L in stack: raise RecursionError
        stack.add(L)
        if not dH[L]: G[L] = base_star if L[0] == 's' else base0
        else:
            top = max(sum(val(y) for y in w) for w in dH[L])
            G[L] = top if L in theta else top + 1
        stack.discard(L); return G[L]
    try:
        for L in H: val(L)
    except RecursionError: return None
    return {L: v for L, v in G.items() if L not in theta}

def part_shape(m, letters, dval, theta, piv, verbose=True):
    H, dH, bad_dg, bad_d2 = h.projection(None, m, letters, dval, theta, piv)
    pi_ok = not bad_dg and not bad_d2
    Th = [L for L in H if L in theta]
    found = None
    for bs, b0 in itertools.product(range(0, 5), range(0, 3)):
        g = g_recursive(m, H, dH, theta, bs, b0)
        if g is None: continue
        r, why = h.shape(m, H, dH, theta, g)
        if r: found = r; break
    return H, dH, pi_ok, Th, found, why

# ---------------- H(I) engine on n vertices (cell loop of workshop/rounds/106/theorist_loopfilt.run, vertices generalised)
P = lf.P
def engine(m, letters, wt, ideal_letters, verts):
    dth = {L: m.partial({(L,): 1}) for L in ideal_letters}
    cells = bad = tot = notclosed = 0; hsum = 0; det = []
    ws_all = {}
    for s in verts:
        for t in verts:
            ws = []
            def rec(cur, v):
                if v == t and cur: ws.append(tuple(cur))
                if len(cur) > 3 * len(verts) + 4: return
                for L in letters:
                    if m.src(L) == v: cur.append(L); rec(cur, m.tgt(L)); cur.pop()
            rec([], s)
            if not ws: continue
            bydeg = {}
            for w in ws: bydeg.setdefault(m.wdeg(w), []).append(w)
            piv_d = {}
            for dg, L in bydeg.items():
                piv = {}
                for w in L:
                    if any(x in dth for x in w): lf.reduce_into(piv, {w: 1})
                piv_d[dg] = piv
            for d1, L in bydeg.items():
                for w in L:
                    for k, x in enumerate(w):
                        if x in dth: lf.reduce_into(piv_d.setdefault(d1 - 1, {}), {w[:k] + u + w[k + 1:]: c for u, c in dth[x].items()})
            rk = {}; dimI = {}
            for dg, piv in piv_d.items():
                dimI[dg] = len(piv); img = {}; r = 0
                for row in piv.values():
                    e = {w: c % P for w, c in m.partial(row).items() if c % P}
                    if e:
                        if lf.reduce_into(dict(piv_d.get(dg - 1, {})), dict(e)) is not None: notclosed += 1
                        if lf.reduce_into(img, e) is not None: r += 1
                rk[dg] = r
            for dg in dimI:
                hh = dimI[dg] - rk.get(dg, 0) - rk.get(dg + 1, 0); cells += 1; tot += dimI[dg]
                if hh: bad += 1; hsum += hh; det.append((s, t, dg, dimI[dg], hh))
    return dict(cells=cells, dimI=tot, notclosed=notclosed, bad_cells=bad, sumH=hsum, first=det[:3])

def lf_line(m, letters, wt):
    verts = sorted({m.src(L) for L in letters} | {m.tgt(L) for L in letters})
    return h.lf_check(verts, [(m.src(L), m.tgt(L), wt(L)) for L in letters])

def h1w(letters, dval, theta):
    """(H1w), the part of (H1) that Lemmas 1-2 use: every Theta letter with a linear term has a pivot x (coefficient != 0,
    x not a Theta letter, x in no OTHER Theta letter's linear part, pivots distinct).  Non-Theta letters may have linear terms."""
    src = [L for L in letters if L in theta and h.lin_part(dval, L)]
    piv = {}; bad = []
    for L in sorted(src, key=repr):
        lp = h.lin_part(dval, L)
        x = next((x for x, k in lp.items() if k and x not in theta and x not in piv.values()
                  and all(x not in h.lin_part(dval, M) for M in src if M != L)), None)
        if x is None: bad.append(h.nm(L))
        else: piv[L] = x
    return bad, piv

def analyse(m, letters, wt, dval, theta, verts, ctrl=False, verbose=False):
    d2 = [h.nm(L) for L in letters if m.partial(dval[L])]
    lfok, mins = lf_line(m, letters, wt)
    cyc = any(v != float('inf') for v in mins.values())
    r1, _ = h.h1(letters, dval, theta)
    bad, piv = h1w(letters, dval, theta)
    out = dict(letters=len(letters), theta=len(theta), d2=len(d2), LF=lfok, cycles=cyc, H1=not any(r1.values()), H1w=not bad,
               nontheta_lin=r1['nontheta_src'], piv=[(h.nm(L), h.nm(x)) for L, x in piv.items()])
    if not bad:
        H, dH, pi_ok, Th, found, why = part_shape(m, letters, dval, theta, piv)
        out.update(PI=pi_ok, thetaprime=[h.nm(L) for L in Th], SH=bool(found) or not Th,
                   lead=[(h.nm(t), h.nm(w[0]) + '·' + h.nm(w[1])) for t, w in (found[1].items() if found else [])], why=None if found else why)
        if verbose:
            for L in H:
                if dH[L]: print(f'     d_H {h.nm(L):10s}{"(Theta)" if L in theta else "       "} = {h.show(dH[L])}')
            if found: print('     g =', {h.nm(k): v for k, v in found[0].items() if v})
        if ctrl and found and Th:   # planted g: lower the largest X value to 0
            G = {k: v for k, v in found[0].items() if k not in theta}; x0 = max(G, key=lambda k: G[k]); G[x0] = 0
            r, why2 = h.shape(m, H, dH, theta, G); out['SHctrl'] = 'NO (caught): ' + why2 if not r else 'YES (NOT caught)'
    e = engine(m, letters, wt, sorted(theta, key=repr), verts)
    out['ENG'] = (e['bad_cells'], e['sumH'], e['notclosed'], e['cells'], e['dimI'])
    if ctrl:
        e1 = engine(m, letters, wt, [L for L in letters if L[0] == 's'], verts)
        out['ENGctrl_stars'] = (e1['bad_cells'], e1['sumH'], e1['notclosed'])
    return out

def step(n, row, op, i, ctrl=False, verbose=False):
    nn, rels, _ = lna(n, row, op)
    arr, d = minimal_model(n, rels, op)
    q, m, letters, wt, dval, theta = kM(arr, d, i)
    return analyse(m, letters, wt, dval, theta, list(range(1, n + 1)), ctrl, verbose), (arr, d)

def tree_example():
    """Not an LNA: y: 1->2, al: 2->3, al2: 2->4, z: 3->5, relations y.al, al.z; vertex 2 has two out-arrows.  Construction 2.2
    model written by hand (p, q of degree 1, overlap c of degree 2); H = A not machine-checked here (d^2 = 0 is)."""
    sys.path.insert(0, 'tests')
    from helpers import path_algebra
    from quivermutation import mutation
    tl = load('tl', 'workshop/rounds/102/toolsmith_litex.py')
    Alg = path_algebra([[1, 2], [2, 3], [2, 4], [3, 5]], [[[1, 2, 3]], [[2, 3, 5]]])
    arr = {'y': (1, 2, 0), 'al': (2, 3, 0), 'al2': (2, 4, 0), 'z': (3, 5, 0), 'p': (1, 3, 1), 'q': (2, 5, 1), 'c': (1, 5, 2)}
    d = {'p': {('y', 'al'): 1}, 'q': {('al', 'z'): 1}, 'c': {('p', 'z'): 1, ('y', 'q'): -1}}
    Q = opp.Q(arr, [], 0); Q.d = d
    q, m, letters, wt, dval, theta = kM(arr, d, 2)
    out = analyse(m, letters, wt, dval, theta, [1, 2, 3, 4, 5], ctrl=True, verbose=True)
    print(f'[TREE] F gate at 2: {mutation.mutationIsPossibleAtVertex(Alg, 2)}, J {tl.perI(Alg, 2) or 0}; parent d^2 bad {[L for L in d if Q.dlin(d[L])]}')
    print(f'[TREE] {out}')

def library_steps(n, row):
    sys.path.insert(0, 'tests')
    from quivermutation import nakayama, mutation, pathAlgebra
    tl = load('tl', 'workshop/rounds/102/toolsmith_litex.py')
    A = nakayama.LinearNakayamaAlgebra(n, row); D = pathAlgebra.dualPathAlgebra(A); res = {}
    for v in range(1, n + 1):
        res[v] = ('F', mutation.mutationIsPossibleAtVertex(A, v), not tl.perI(A, v), 'R', mutation.mutationIsPossibleAtVertex(D, v), not tl.perI(D, v))
    return res

FEATURED = [(5, (2, 2, 0), 1), (6, (3, 3, 0, 0), 1)]   # Def. 3.2 on A at i = the repo's F step at i (E-252: F = Oppermann's left)
ROWS = [(4, (2, 0)), (4, (3, 0)), (5, (2, 2, 0)), (5, (2, 3, 0)), (5, (3, 2, 0)), (5, (0, 2, 0)), (5, (2, 0, 2)), (6, (2, 2, 2, 0)), (6, (3, 3, 0, 0)), (6, (2, 4, 0, 0))]

if __name__ == '__main__':
    rows = ROWS
    if '--rows' in sys.argv:
        rows = [(int(x.split(':')[0]), tuple(int(c) for c in x.split(':')[1])) for x in sys.argv[sys.argv.index('--rows') + 1].split(',')]
    # [LF] positive parent weights alone do not give (LF): kQ/(a1 c), a1: i->1, c: 1->2, a2: i->2, f: i->2 (deg 1, d f = a1 c)
    for w2 in (5, 1):
        arr = {'a1': (0, 1, 0), 'c': (1, 2, 0), 'a2': (0, 2, 0), 'f': (0, 2, 1)}
        q = lf.mk(arr, ['a1', 'a2'], {'f': {('a1', 'c'): 1}}, {'a1': 1, 'c': 1, 'a2': w2, 'f': 2})
        m = opp.M(q, {'last': -1}); letters = opp.all_letters(q) + [('s', al) for al in sorted(q.A)]   # lf.setup would assert wt >= 1
        wt = lambda L: -q.w[L[1]] if L[0] == 's' else q.w[L[2]] + (q.w[L[1]] if L[1] else 0) - (q.w[L[3]] if L[3] else 0)
        ok, mins = lf_line(m, letters, wt)
        print(f'[LF] two A-arrows, w(a1) = 1, w(a2) = {w2}, w(f) = 2: min cycle weight per vertex {mins} -> {"YES" if ok else "NO"}'
              f' (condition (C) w(f) > max w(A): {2 > max(1, w2)})')
    ok, mins = h.lf_check([0, 1], [(0, 1, 2), (1, 0, -2)])
    print(f'[LF control] toy quiver with a weight-0 cycle -> {"YES (NOT caught)" if ok else "NO (caught)"}')
    # [MM] controls
    nn, rels, op = lna(5, (2, 2, 0))
    arr, d = minimal_model(nn, rels, op, plant=True)
    print('[MM control] A5 22 0 with one term of one d dropped: failures', len(check_mm(nn, rels, op, arr, d)), '-> must be > 0')
    for n, row in rows:
        lib = library_steps(n, row)
        for op in (False, True):
            nn, rels, _ = lna(n, row, op)
            arr, d = minimal_model(n, rels, op)
            mm = check_mm(n, rels, op, arr, d)
            print(f'=== A{n} row {"".join(map(str, row))} {"op" if op else "  "}: parent letters {len(arr)} (deg>=1: {sum(1 for v in arr.values() if v[2])}), [MM] H = A: {"YES" if not mm else "NO " + str(mm[:3])}')
            for i in range(1, n + 1):
                if not any(s == i and g == 0 for s, t, g in arr.values()): continue   # sink: no A arrow, kM = parent
                feat = not op and (n, row, i) in FEATURED
                out, _ = step(n, row, op, i, ctrl=feat, verbose=feat)
                F, R = lib[i][1:3], lib[i][4:6]
                tag = f'F gate {int(F[0])} J0 {int(F[1])} | R gate {int(R[0])} J0 {int(R[1])}'
                print(f'  i={i} [{tag}] letters {out["letters"]} Theta {out["theta"]} D2bad {out["d2"]} LF {out["LF"]} (cycles {out["cycles"]}) '
                      f'H1 {out["H1"]} (non-Theta linear {out["nontheta_lin"]}) H1w {out["H1w"]} piv {out["piv"]} PI {out.get("PI")} Theta\' {out.get("thetaprime")} SH {out.get("SH")} '
                      f'lead {out.get("lead")} ENG(bad,sumH,notclosed,cells,dimI) {out["ENG"]}' + (f' why {out["why"]}' if out.get('why') and not out.get('SH') else ''))
                for k in ('SHctrl', 'ENGctrl_stars'):
                    if k in out: print(f'     control {k}: {out[k]}')
    tree_example()
