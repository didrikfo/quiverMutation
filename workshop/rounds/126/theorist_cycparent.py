"""Round 126 (theorist), T5 / H-023, H-031: the hypotheses (LF), (H1w), (H3) and H(I) = 0 on REAL Construction-2.2 parents
with a loop letter at the mutation vertex: cyclic Nakayama algebras kZ_n / (zero relations), path-length weights.
Parent = minimal model built Tate-style cell by cell, cells (source s, weight W) with target s + W mod n (exact rationals);
H = A is re-checked by an independent pass.  kM = Def. 3.2 at i (workshop/rounds/069, signs last = -1), checks from
workshop/rounds/123/theorist_lnaparent.py and workshop/rounds/122/theorist_h031.py.  H(I) is computed here per cell
(s, t, shifted weight, degree) with the (C)-shift phi(i) = c, every letter weight > 0, so every cell is finite: no word cap.
Also tests the candidate g = "path length" for (H3): g_lam(x) = 2 * shifted weight + lam * degree, lam = 1..3.
Controls: [MM] dropped term must fail H = A; [ENG] closed ideal (stars) must give H != 0; [SH] planted g must fail;
[LF] the 2-cycle with w(b) = -1 shift in place of the (C) shift must show a nonpositive letter.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/126/theorist_cycparent.py [WMAX]"""
import sys, importlib.util, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '.')
def load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
P = load('lnap', 'workshop/rounds/123/theorist_lnaparent.py'); h = P.h; lf = P.lf; opp = P.opp

def nonzero(n, rels, s, W):
    """path of length W from s on the cycle 1..n is nonzero iff it contains no relation path (a, l)"""
    for a, l in rels:
        for o in range(0, W - l + 1):
            if (s - 1 + o) % n == (a - 1) % n: return False
    return True

def kupisch(n, rels):
    """Kupisch series: c[j] = length of P_(j+1) = min over relations (s, l) of ((s - 1 - j) mod n + l)"""
    return tuple(min((s - 1 - j) % n + l for s, l in rels) for j in range(n))

def gldim_finite(c):
    """gl.dim < inf iff every rad P_j has a terminating uniserial syzygy chain: Omega(top t, length l) = (t + l, c_t - l)"""
    n = len(c)
    for j in range(n):
        t, l, seen = (j + 1) % n, c[j] - 1, set()
        while 0 < l < c[t]:
            if (t, l) in seen: return False
            seen.add((t, l)); t, l = (t + l) % n, c[t] - l
    return True

def kupisch_classes(n, LMAX):
    """independent enumeration: all admissible Kupisch series (c_(j+1) >= c_j - 1, c_j >= 2) whose relations (starts j with
    c_(j+1) >= c_j, length c_j) have length <= LMAX, finite gl.dim, one per rotation"""
    out = set()
    for c in itertools.product(range(2, n + LMAX), repeat=n):
        if any(c[(j + 1) % n] < c[j] - 1 for j in range(n)): continue
        if max(c[j] for j in range(n) if c[(j + 1) % n] >= c[j]) > LMAX or not gldim_finite(c): continue
        out.add(min(c[k:] + c[:k] for k in range(n)))
    return out

def minimal_model(n, rels, WMAX, plant=False, cap=12):
    """stops early once LMAX consecutive weights add no letter (complete: an Anick (m+1)-chain extends an m-chain by fewer
    than LMAX arrows), or once more than cap letters of degree >= 1 exist (reported as incomplete)"""
    arr = {}; d = {}; wt = {}; LM = max(l for _, l in rels)
    tg = lambda s, W: (s - 1 + W) % n + 1
    q = opp.Q(arr, [], 0); q.d = d
    def words(s, W):
        out = []
        def rec(v, w, cur):
            if w == W and len(cur) >= 2: out.append(tuple(cur)); return
            for L, (a, b, _) in arr.items():
                if a == v and w + wt[L] <= W and wt[L] < W: rec(b, w + wt[L], cur + [L])
        rec(s, 0, []); return out
    cnt = 0
    for W in range(1, WMAX + 1):
        for s in range(1, n + 1):
            t = tg(s, W)
            if W == 1: arr['x%d' % s] = (s, t, 0); wt['x%d' % s] = 1; continue
            bydeg = {}
            for w in words(s, W): bydeg.setdefault(q.deg(w), []).append(w)
            for m in sorted(bydeg):
                Z = P.kernel(bydeg[m], [q.dpath(w) for w in bydeg[m]]) if m > 0 else [{w: 1} for w in bydeg[m]]
                B = {}
                for w in bydeg.get(m + 1, []): P.rref_add(B, q.dpath(w))
                reps = [z for z in Z if P.rref_add(B, z)]
                if m == 0 and nonzero(n, rels, s, W): reps = reps[1:]
                for z in reps:
                    cnt += 1; L = ('r%d' if m == 0 else 'c%d_%%d' % (m + 1)) % cnt
                    arr[L] = (s, t, m + 1); d[L] = P.integral(z); wt[L] = W
        if W >= max(wt.values()) + LM or cnt > cap: break
    if plant:
        L = next(L for L in d if len(d[L]) >= 2); d[L].pop(next(iter(d[L])))
    return arr, d, wt

def check_mm(n, rels, arr, d, wt, WMAX):
    q = opp.Q(arr, [], 0); q.d = d
    bad = [L for L in d if q.dlin(d[L])]
    for W in range(1, WMAX + 1):
        for s in range(1, n + 1):
            ws = []
            def rec(v, w, cur):
                if w == W and cur: ws.append(tuple(cur)); return
                for L, (a, b, _) in arr.items():
                    if a == v and w + wt[L] <= W: rec(b, w + wt[L], cur + [L])
            rec(s, 0, [])
            bydeg = {}
            for w in ws: bydeg.setdefault(q.deg(w), []).append(w)
            for m in sorted(set(bydeg) | {0}):
                z = len(P.kernel(bydeg.get(m, []), [q.dpath(w) for w in bydeg.get(m, [])])) if m > 0 else len(bydeg.get(0, []))
                Bb = {}
                for w in bydeg.get(m + 1, []): P.rref_add(Bb, q.dpath(w))
                want = 1 if (m == 0 and nonzero(n, rels, s, W)) else 0
                if z - len(Bb) != want: bad.append((s, W, m, z - len(Bb), want))
    return bad

def kM(arr, d, wt, i):
    A = sorted(L for L, (s, t, g) in arr.items() if s == i and g == 0)
    q = lf.mk(dict(arr), A, {k: dict(v) for k, v in d.items()}, dict(wt)); q.i = i
    m, letters, wtf = lf.setup(q)
    dval = {L: {w: k for w, k in m.partial({(L,): 1}).items() if k} for L in letters}
    theta = {L for L in letters if L[0] == 'm' and L[3] is not None}
    return q, m, letters, wtf, dval, theta, A

def shifted(arr, wt, i, m, letters, wtf, c=None):
    A = [L for L, (s, t, g) in arr.items() if s == i and g == 0]
    nonA = [wt[L] for L, (s, t, g) in arr.items() if s == i and L not in A]
    mA = max(wt[a] for a in A); mN = min(nonA, default=mA + 2)
    C = mA < mN
    if c is None: c = Fr(mA + mN, 2)
    phi = lambda v: c if v == i else 0
    return C, {L: 2 * (wtf(L) + phi(m.tgt(L)) - phi(m.src(L))) for L in letters}   # doubled: integers

def engine(m, letters, sw, ideal, verts, WMAX2):
    """H(I) per cell (s, t, shifted weight <= WMAX2, degree), I = (ideal, d ideal); exact mod p"""
    dth = {L: m.partial({(L,): 1}) for L in ideal}
    res = dict(cells=0, dimI=0, notclosed=0, bad=0, sumH=0, first=[])
    for s in verts:
        cells = {}
        def rec(cur, v, w):
            if cur: cells.setdefault((v, w), []).append(tuple(cur))
            for L in letters:
                if m.src(L) == v and w + sw[L] <= WMAX2: cur.append(L); rec(cur, m.tgt(L), w + sw[L]); cur.pop()
        rec([], s, 0)
        for (t, W), ws in cells.items():
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
                        if x in dth: lf.reduce_into(piv_d.setdefault(d1 - 1, {}), {w[:k] + u + w[k + 1:]: cc for u, cc in dth[x].items()})
            rk = {}; dimI = {}
            for dg, piv in piv_d.items():
                dimI[dg] = len(piv); img = {}; r = 0
                for row in piv.values():
                    e = {w: cc % lf.P for w, cc in m.partial(row).items() if cc % lf.P}
                    if e:
                        if lf.reduce_into(dict(piv_d.get(dg - 1, {})), dict(e)) is not None: res['notclosed'] += 1
                        if lf.reduce_into(img, e) is not None: r += 1
                rk[dg] = r
            for dg in dimI:
                hh = dimI[dg] - rk.get(dg, 0) - rk.get(dg + 1, 0); res['cells'] += 1; res['dimI'] += dimI[dg]
                if hh: res['bad'] += 1; res['sumH'] += hh; res['first'] = (res['first'] + [(s, t, W, dg, hh)])[:3]
    return res

def g_pathlength(m, H, dH, theta, sw, lam):
    """candidate g = shifted weight (doubled) + lam * degree on X; Theta values are set by shape()"""
    return {L: sw[L] + lam * m.deg(L) for L in H if L not in theta}

def analyse(n, rels, arr, d, wt, i, WENG, ctrl=False, verbose=False):
    q, m, letters, wtf, dval, theta, A = kM(arr, d, wt, i)
    out = dict(letters=len(letters), theta=len(theta))
    out['d2'] = sum(1 for L in letters if m.partial(dval[L]))
    out['nonhom'] = sum(1 for L in letters for w in dval[L] if sum(wtf(x) for x in w) != wtf(L))
    loops = [h.nm(L) for L in letters if L[0] == 'm' and m.src(L) == m.tgt(L)]
    out['parent_loops_at_i'] = [L for L, (s, t, g) in arr.items() if s == i and t == i]
    lfok, mins = P.lf_line(m, letters, wtf); out['LF'] = lfok; out['cycles'] = any(v != float('inf') for v in mins.values())
    out['kM_loop_letters'] = len(loops)
    C, sw = shifted(arr, wt, i, m, letters, wtf); out['C'] = C; out['minshift'] = Fr(min(sw.values()), 2)
    bad, piv = P.h1w(letters, dval, theta); out['H1w'] = not bad
    r1, _ = h.h1(letters, dval, theta); out['H1strict'] = not any(r1.values())
    if not bad:
        H, dH, pi_ok, Th, found, why = P.part_shape(m, letters, dval, theta, piv)
        out['PI'] = pi_ok; out['thetaprime'] = [h.nm(L) for L in Th]; out['SH'] = bool(found) or not Th
        out['lead'] = [(h.nm(t), h.nm(w[0]) + '·' + h.nm(w[1])) for t, w in (found[1].items() if found else [])]
        if not found and Th: out['why'] = why
        pl = {}
        for lam in (1, 2, 3):
            r, why2 = h.shape(m, H, dH, theta, g_pathlength(m, H, dH, theta, sw, lam)); pl[lam] = 'YES' if r else why2
        out['g=pathlength'] = pl if Th else 'vacuous'
        if verbose:
            for L in H:
                if dH[L]: print(f'     d_H {h.nm(L):12s}{"(Theta)" if L in theta else "       "} = {h.show(dH[L])}')
            if found: print('     g (recursive) =', {h.nm(k): v for k, v in found[0].items() if v})
        if ctrl and found and Th:
            G = {k: v for k, v in found[0].items() if k not in theta}; x0 = max(G, key=lambda k: G[k]); G[x0] = 0
            r, why2 = h.shape(m, H, dH, theta, G); out['SHctrl'] = 'NO (caught): ' + why2 if not r else 'YES (NOT caught)'
    verts = list(range(1, n + 1))
    e = engine(m, letters, sw, sorted(theta, key=repr), verts, WENG)
    out['ENG'] = (e['bad'], e['sumH'], e['notclosed'], e['cells'], e['dimI'])
    if ctrl:
        e1 = engine(m, letters, sw, [L for L in letters if L[0] == 's'], verts, WENG)
        out['ENGctrl_stars'] = (e1['bad'], e1['sumH'], e1['first'][:1])
    return out

# cyclic Nakayama parents: (n, relations (start, length)); minimal sets of zero relations, finite global dimension
ROWS = [(2, [(1, 2)]), (2, [(1, 3)]), (3, [(1, 2)]), (3, [(1, 3)]), (3, [(1, 2), (2, 2)]), (3, [(1, 4)]), (2, [(1, 2), (2, 4)]),
        (3, [(1, 3), (2, 3)]), (4, [(1, 2), (2, 2)]), (4, [(1, 3), (3, 3)]), (4, [(1, 4), (2, 4)]), (3, [(1, 2), (2, 3)])]

if __name__ == '__main__':
    WMAX = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 14
    sel = ROWS
    if '--sweep' in sys.argv:   # all antichains of zero relations (s, l), one per start, 2 <= l <= LMAX, on Z_n, finite gl.dim only
        nn, LMAX = (int(x) for x in sys.argv[sys.argv.index('--sweep') + 1].split(':'))
        sel = []
        for ls in itertools.product([0] + list(range(2, LMAX + 1)), repeat=nn):
            rels = [(s0 + 1, l) for s0, l in enumerate(ls) if l]
            if not rels or ls != max(ls[k:] + ls[:k] for k in range(nn)): continue   # one representative per rotation
            # antichain: no relation path contains another
            def inside(r1, r2):   # r1 inside r2 ?
                return not nonzero(nn, [r1], r2[0], r2[1]) and r1 != r2
            if any(inside(a, b) for a in rels for b in rels): continue
            sel.append((nn, rels))
    CAP = int(sys.argv[sys.argv.index('--cap') + 1]) if '--cap' in sys.argv else 60
    if '--rows' in sys.argv:
        sel = [ROWS[int(k)] for k in sys.argv[sys.argv.index('--rows') + 1].split(',')]
    arr, d, wt = minimal_model(2, [(1, 2), (2, 2)], 10)
    print(f'[GLDIM control] 2-cycle x1x2 = x2x1 = 0 (selfinjective, infinite gl.dim): letters up to W = 10: {len(arr)}, max weight {max(wt.values())} (must grow with W); Kupisch gl.dim finite: {gldim_finite(kupisch(2, [(1, 2), (2, 2)]))} (must be False)')
    arr, d, wt = minimal_model(3, [(1, 2), (2, 2)], WMAX, plant=True)
    print(f'[MM control] Z3/(x1x2, x2x3) with one dropped term: H = A failures {len(check_mm(3, [(1, 2), (2, 2)], arr, d, wt, WMAX))} (must be > 0)')
    tally = dict(classes=0, inf=0, inf_but_complete=0, run=0, capskip=0, inst=0)
    runk = set()
    for k, (n, rels) in enumerate(sel):
        tally['classes'] += 1
        kc = kupisch(n, rels)
        if not gldim_finite(kc):   # skip on a PROVED property; cross-check: the cap-12 builder must not find a complete model
            arr, d, wt = minimal_model(n, rels, WMAX, cap=12); top = max(wt.values())
            comp = top <= WMAX - max(l for _, l in rels) and sum(1 for v in arr.values() if v[2]) <= 12
            tally['inf'] += 1; tally['inf_but_complete'] += comp
            print(f'=== Z{n} rels {rels}: SKIPPED, infinite gl.dim (Kupisch {kc}); cap-12 builder complete: {comp} (must be False)'); continue
        arr, d, wt = minimal_model(n, rels, WMAX, cap=CAP)
        top = max(wt.values())
        fin = top <= WMAX - max(l for _, l in rels) and sum(1 for v in arr.values() if v[2]) <= CAP
        if not fin:
            tally['capskip'] += 1
            print(f'=== Z{n} rels {rels}: SKIPPED BY SIZE CAP (finite gl.dim, Kupisch {kc}): model not complete (more than {CAP} letters of degree >= 1 or weight > WMAX - LMAX)'); continue
        tally['run'] += 1; tally['inst'] += n; runk.add(min(kc[j:] + kc[:j] for j in range(n)))
        mm = check_mm(n, rels, arr, d, wt, min(WMAX, top + 2 * n + 2))
        print(f'=== Z{n} rels {rels}: letters {len(arr)}, degrees {sorted(g for *_, g in arr.values())}, max weight {top} '
              f'(complete: {fin}), [MM] H = A to W = {min(WMAX, top + 2 * n + 2)}: {"YES" if not mm else "NO " + str(mm[:3])}')
        print('    parent d:', {L: h.show({tuple(w): c for w, c in d[L].items()}) for L in d})
        for i in range(1, n + 1):
            feat = (n, rels) == (3, [(1, 2), (2, 2)]) and i == 1
            out = analyse(n, rels, arr, d, wt, i, WENG=2 * (top + 6), ctrl=feat, verbose=feat)
            print(f'  i={i} ' + ' '.join(f'{k} {v}' for k, v in out.items()))
    print('[TALLY]', tally)
    if '--sweep' in sys.argv:
        kk = kupisch_classes(nn, LMAX)
        print(f'[KUPISCH control] finite-gl.dim classes for {nn}:{LMAX} by Kupisch series: {len(kk)}; run: {len(runk)}; '
              f'missing from run: {sorted(kk - runk)}; run but not in list: {sorted(runk - kk)}')
    # [LF] control: the M1 parent with the shift c = max w(A) (outside the (C) interval) must give a nonpositive letter
    arr, d, wt = minimal_model(2, [(1, 2)], 8)
    q, m, letters, wtf, dval, theta, A = kM(arr, d, wt, 1)
    _, sw = shifted(arr, wt, 1, m, letters, wtf, c=1)
    print(f'[LF control] Z2/(x1x2) at i = 1, shift c = max w(A) = 1: nonpositive letters {[h.nm(L) for L in letters if sw[L] <= 0]} (must be non-empty)')
