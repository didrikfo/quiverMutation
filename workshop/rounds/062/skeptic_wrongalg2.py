"""Round 062 (skeptic), T10 / E-177 open item: equal-dimension wrong-algebra control for the E-180 compare2 dimension assertion at
 (a) step 12 of class 1 (child with TWO parallel pairs P1=(5,7), P2=(7,4)) and (b) the 9 class-2 failing steps (no parallel arrows).
Needs /tmp/tsm/c1.pkl, /tmp/tsm/c2.pkl; rebuild each (twice, resumable; ~8 and ~6 min):
  DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 CLS 20000 /tmp/tsm/cCLS.pkl 100
Usage (repo root):  CLS=1 timeout 10m .venv/bin/python workshop/rounds/062/skeptic_wrongalg2.py     (step 12; MAXV env caps variants)
                    CLS=2 timeout 10m .venv/bin/python workshop/rounds/062/skeptic_wrongalg2.py     (9 steps)
Wrong algebras = child quiver with a changed relation set; kept only if dim K Q/I' = d pair by pair (own linear algebra, quotient_dims).
Verdict = compare-style: dimcheck (E-180) then symcheck2 (E-174).  TRUTH is independent of both:
 class 1 step 12: relations are 3 'line' relations x*(l0 b4 + l1 b5)=0 for x in {1->5, 3->5, 6->5} (points p1,p3,p7 of P(P1)) and the kernel of a bilinear
   form phi on P1 (x) P2 (3-dim, = the 3 relations on 5->7->4).  Group GL2(P1) x GL2(P2) x scalars; no arrow has a longer parallel path, so this is all of
   Aut over the quiver.  Orbits: phi rank 2: coincidence pattern of (p1,p3,p7) [GL2(P2) transitive on nondegenerate phi]; phi = alpha (x) beta rank 1:
   pattern of (p1,p3,p7,q) with q = ker alpha.  Iso to true child iff same orbit (true: p = (0:1),(1:0),(1:1) distinct, phi rank 2).
 class 2: no arrow is parallel to a longer path, so Aut is the torus (arrow -> t_a arrow); truth = exists t in mu_4^arrows with t(I) = I' (checked in GF(13),
   4th roots 1,5,12,8; first tried t in {+-1}), comparing the ideals as subspaces of every path space."""
import pickle, itertools, collections, sys, os, time
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
from fractions import Fraction
MAXV = int(os.environ.get('MAXV', '100000')); P = 13

def paths_of(arrows):
    out = collections.defaultdict(list)
    for a in arrows: out[a[0]].append(a)
    verts = sorted({a[0] for a in arrows} | {a[1] for a in arrows})
    paths = collections.defaultdict(list)
    for v in verts: paths[(v, v)].append(())
    fr = [((a,), a[0], a[1]) for a in arrows]
    while fr:
        nx = []
        for p, s, t in fr:
            paths[(s, t)].append(p)
            for a in out[t]: nx.append((p + (a,), s, a[1]))
        fr = nx
    return paths

def ideal_space(arrows, rels, paths, tv=None):
    """per (s,t): RREF (mod P) of the span of pre*r*post, with r scaled by the torus element tv (dict arrow->int mod P)."""
    res = {}
    for (s, t), ps in paths.items():
        idx = {p: i for i, p in enumerate(ps)}; rows = []
        for r in rels:
            r_s = next(iter(r))[0][0]; r_t = next(iter(r))[-1][1]
            for pre in paths.get((s, r_s), []):
                for post in paths.get((r_t, t), []):
                    row = [0] * len(ps)
                    for rp, cf in r.items():
                        q = pre + tuple(rp) + post
                        if q in idx:
                            sc = 1
                            if tv:
                                for a in rp: sc = sc * tv[a] % P
                            cf = Fraction(cf); row[idx[q]] = (row[idx[q]] + cf.numerator * pow(cf.denominator, -1, P) * sc) % P
                    if any(row): rows.append(row)
        # RREF
        rk = 0; M = rows
        for c in range(len(ps)):
            pr = next((i for i in range(rk, len(M)) if M[i][c]), None)
            if pr is None: continue
            M[rk], M[pr] = M[pr], M[rk]; inv = pow(M[rk][c], -1, P); M[rk] = [x * inv % P for x in M[rk]]
            for i in range(len(M)):
                if i != rk and M[i][c]: f = M[i][c]; M[i] = [(x - f * y) % P for x, y in zip(M[i], M[rk])]
            rk += 1
        res[(s, t)] = tuple(tuple(r) for r in M[:rk])
    return res

def torus_iso(arrows, true_rels, var_rels):
    paths = paths_of(arrows); target = ideal_space(arrows, var_rels, paths)
    used = sorted({a for r in true_rels for p in r for a in p} | {a for r in var_rels for p in r for a in p})
    for vals in ([1, P - 1], [1, 5, P - 1, 8]):
        for ch in itertools.product(vals, repeat=len(used)):
            tv = dict(zip(used, ch))
            if ideal_space(arrows, true_rels, paths, tv) == target: return True
    return False

def no_parallel_longer(arrows):
    paths = paths_of(arrows)
    return all(len([p for p in paths[(a[0], a[1])] if len(p) >= 2]) == 0 for a in arrows)

def no_parallel_arrows(arrows):
    return len({(a[0], a[1]) for a in arrows}) == len(arrows) and all(a[0] != a[1] for a in arrows)

def torus_iso_full(arrows, true_rels, var_rels):
    """Full torus (not only mu_4). The true ideal = monomials + exactly ONE binomial c1*p + c2*q; the torus image of it is monomials +
    c1*p + c2*lam*q with lam = t^q/t^p; every lam in GF(13)* is tried (a superset of what the torus reaches), so NOT here is a proof of NOT."""
    bins = [r for r in true_rels if len(r) == 2]; assert len(bins) == 1, 'need exactly one binomial'
    paths = paths_of(arrows); target = ideal_space(arrows, var_rels, paths); (p, c1), (q, c2) = list(bins[0].items())
    for lam in range(1, P):
        rels = [({p: c1, q: Fraction(c2) * lam} if r is bins[0] else r) for r in true_rels]
        if ideal_space(arrows, rels, paths) == target: return True
    return False

def tally(verd, eqdim, truth, vd):
    return ('eqdim' if eqdim else 'DIMDIFF', 'truth=iso' if truth else 'truth=NOT', 'test=' + ('iso' if vd.startswith('iso') else 'NO'))

recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
mu4 = collections.Counter(); byrank = collections.Counter(); table = collections.Counter(); strings = collections.Counter(); bad = []

def run_variant(n_, TE, arrows, carrows, crels, d, new, truth, tag):
    qd = quotient_dims(carrows, new)
    eqdim = all(qd.get(k, 0) == d.get(k, 0) for k in set(qd) | set(d) if k[0] != k[1])
    if not eqdim: table[('DIMDIFF', 'skipped')] += 1; return
    dc = dimcheck(carrows, new, d)
    vd = symcheck2(TE, arrows, 'same', carrows, new)
    full = 'NO' if dc else ('iso' if vd.startswith('iso') else 'NO')
    table[('eqdim', 'truth=iso' if truth else 'truth=NOT', 'dimcheck=' + ('rej' if dc else 'pass'), 'symcheck2=' + ('iso' if vd.startswith('iso') else 'NO'))] += 1
    strings[('truth=iso' if truth else 'truth=NOT', vd.split(' (')[0] if vd.startswith('NO') else 'iso')] += 1
    if (full == 'iso') != truth: bad.append((n_, tag, truth, vd))

if CLS == 1:
    n_ = 12; r = recs[n_]; a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); A, d, carrows, crels = child_info(c)
    assert no_parallel_longer(carrows), 'arrow parallel to longer path'
    b4, b5, c0, c1 = (5, 7, 4), (5, 7, 5), (7, 4, 0), (7, 4, 1)
    L = [i for i, x in enumerate(crels) if all(len(p) == 2 and p[1] in (b4, b5) for p in x)]; assert len(L) == 3, L
    M_ = [i for i, x in enumerate(crels) if all(p[0] in (b4, b5) and p[1] in (c0, c1) for p in x)]; assert len(M_) == 3, M_
    xs = [next(iter(crels[i]))[0] for i in L]; print('step 12 line rels', L, 'x', xs, 'bilinear rels', M_)
    LINES = [(1, 0), (0, 1), (1, 1), (1, -1)]
    def proj_eq(u, w): return u[0] * w[1] - u[1] * w[0] == 0
    def pat(ps): return tuple(next(j for j in range(len(ps)) if proj_eq(ps[j], ps[i])) for i in range(len(ps)))
    def ptof(rl): x = next(iter(rl))[0]; return (rl.get((x, b4), 0), rl.get((x, b5), 0))
    true_pts = [ptof(crels[i]) for i in L]
    # true phi: kernel of true relations; recover F by nullspace of the rel vectors on basis b_k c_l
    basis = [(b4, c0), (b4, c1), (b5, c0), (b5, c1)]
    def kernel_rels(F):                      # F[k][l] on b_k c_l (k: 0 -> b4, 1 -> b5; l: 0 -> c0, 1 -> c1); relations = basis of ker phi
        row = sympy.Matrix([[F[0][0], F[0][1], F[1][0], F[1][1]]]); out = []
        for ns in row.nullspace():
            den = sympy.ilcm(*[sympy.Rational(x).q for x in ns]); out.append({basis[j]: Fraction(int(ns[j] * den)) for j in range(4) if ns[j] != 0})
        return out
    # true phi: nullspace of the 3 true bilinear relations
    Mt = sympy.Matrix([[sympy.Rational(int(crels[i].get(bs, 0))) for bs in basis] for i in M_]); phi_true = Mt.nullspace(); assert len(phi_true) == 1
    Ft = [[phi_true[0][0], phi_true[0][1]], [phi_true[0][2], phi_true[0][3]]]; print('true phi', Ft, 'true points', true_pts)
    assert pat(true_pts) == (0, 1, 2) and sympy.Matrix(Ft).det() != 0
    def rank_phi(F): return sympy.Matrix(F).rank()
    def truth_of(pts, F):
        rk = rank_phi(F)
        if rk == 2: return pat(pts) == pat(true_pts)
        # rank 1: F = alpha (x) beta: alpha = nonzero row-combination -> alpha_k = F[k][l] for some l with a nonzero column
        col = next(l for l in range(2) if F[0][l] != 0 or F[1][l] != 0); al = (F[0][col], F[1][col]); q = (al[1], -al[0])
        return False       # rank 1 cannot be iso to a rank 2 phi (invariant: rank of phi on P1 (x) P2)
    Fs = []
    for e in itertools.product([-1, 0, 1], repeat=4):
        if any(e):
            f = tuple(e); neg = tuple(-x for x in e)
            if neg in Fs: continue
            Fs.append(f)
    nv = 0; t0 = time.time()
    for F in Fs:
        Fm = [[F[0], F[1]], [F[2], F[3]]]; kr = kernel_rels(Fm)
        for choice in itertools.product(range(4), repeat=3):
            if nv >= MAXV: break
            new = list(crels); pts = [LINES[ch] for ch in choice]
            for i, x, p_ in zip(L, xs, pts): new[i] = {t: Fraction(w) for t, w in {(x, b4): p_[0], (x, b5): p_[1]}.items() if w != 0}
            for i, kr_ in zip(M_, kr): new[i] = kr_
            if len(kr) != 3: continue
            nv += 1
            n0 = sum(table.values()) - table[('DIMDIFF', 'skipped')]
            run_variant(12, TE, arrows, carrows, crels, d, new, truth_of(pts, Fm), (F, choice))
            if sum(table.values()) - table[('DIMDIFF', 'skipped')] > n0: byrank[(rank_phi(Fm), 'pattern ' + str(pat(pts)))] += 1
        if nv >= MAXV: break
        print('phi', F, 'rank', rank_phi(Fm), 'done, variants so far', nv, 'table', sum(table.values()), '%.0fs' % (time.time() - t0), flush=True)
else:
    for n_, r in enumerate(recs):
        a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
        TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); A, d, carrows, crels = child_info(c)
        assert no_parallel_longer(carrows), 'arrow parallel to longer path'
        assert no_parallel_arrows(carrows), 'two arrows share endpoints'   # needed for Aut = torus
        paths = paths_of(carrows); pool = []
        for (s, t), ps in paths.items():
            long = [p for p in ps if len(p) >= 2]
            for p in long: pool.append({p: Fraction(1)})
            for p, q in itertools.combinations(long, 2):
                pool.append({p: Fraction(1), q: Fraction(-1)})
        seen2 = set(); cnt = collections.Counter(); t0 = time.time()
        for i in range(len(crels)):
            for rep in pool:
                new = crels[:i] + [rep] + crels[i + 1:]
                key = (i, tuple(sorted(tuple(sorted(x.items())) for x in new)))
                if rep == crels[i]: continue
                new_s = ideal_space(carrows, new, paths); tag = (i, tuple(rep.items()))
                qd = quotient_dims(carrows, new)
                if not all(qd.get(k, 0) == d.get(k, 0) for k in set(qd) | set(d) if k[0] != k[1]): table[('DIMDIFF', 'skipped')] += 1; continue
                kk = tuple(sorted(new_s.items()))
                if kk in seen2: continue
                seen2.add(kk)
                truth = torus_iso_full(carrows, crels, new); mu4[(truth, torus_iso(carrows, crels, new))] += 1
                run_variant(n_, TE, arrows, carrows, crels, d, new, truth, tag); cnt[truth] += 1
        # two relations replaced at once (i < j), from the same pool; deduplicated by the resulting ideal
        if os.environ.get('PAIRS', '1') == '1':
            for i, j in itertools.combinations(range(len(crels)), 2):
                for rep1 in pool:
                    for rep2 in pool:
                        if rep1 == crels[i] and rep2 == crels[j]: continue
                        new = list(crels); new[i] = rep1; new[j] = rep2
                        qd = quotient_dims(carrows, new)
                        if not all(qd.get(k, 0) == d.get(k, 0) for k in set(qd) | set(d) if k[0] != k[1]): table[('DIMDIFF', 'skipped')] += 1; continue
                        sp = ideal_space(carrows, new, paths); kk = tuple(sorted(sp.items()))
                        if kk in seen2: continue
                        seen2.add(kk)
                        truth = torus_iso_full(carrows, crels, new); mu4[(truth, torus_iso(carrows, crels, new))] += 1
                        run_variant(n_, TE, arrows, carrows, crels, d, new, truth, ('pair', i, j)); cnt[('p', truth)] += 1
        # sign-flipped binomials (iso controls)
        for i, x in enumerate(crels):
            if len(x) == 2:
                ks = list(x); new = crels[:i] + [{ks[0]: x[ks[0]], ks[1]: -x[ks[1]]}] + crels[i + 1:]
                kk = tuple(sorted(ideal_space(carrows, new, paths).items()))
                if kk in seen2: continue
                seen2.add(kk); tf = torus_iso_full(carrows, crels, new); mu4[(tf, torus_iso(carrows, crels, new))] += 1
                run_variant(n_, TE, arrows, carrows, crels, d, new, tf, ('flip', i)); cnt['flip'] += 1
        print('step', n_, 'v', v, 'eq-dim variants truth iso/NOT', cnt[True], cnt[False], 'flips', cnt['flip'], 'pairs iso/NOT', cnt[('p', True)], cnt[('p', False)], '%.0fs' % (time.time() - t0), flush=True)
print('CLS', CLS)
for k, x in sorted(table.items(), key=str): print(k, x)
print('equal-dim variants by (rank phi, point pattern):', dict(sorted(byrank.items(), key=str)))
print('verdict strings (truth, test):', dict(strings))
print('(full-torus truth, mu_4 truth) counts over distinct equal-dim ideals:', dict(mu4), ' [first ones dedup: singles, pairs, flips share one seen-set per step]')
print('DISAGREEMENTS (full verdict vs truth):', len(bad)); [print('  ', b) for b in bad[:20]]
