"""Round 117 (experimentalist): H-026 K-threshold, even n, power control, class level.  Repo root: .venv/bin/python workshop/rounds/117/experimentalist_kpower.py MODE ...
 table N K1 H1,H2,..   K-table (K = free head stretch, image keys = Coxeter keys of the vertex-1 deletion) for the key classes of lone-3 rows h (one pass, own enumerator
                       copied from rounds/115/experimentalist_kloc.py; libdel optional cross-check against pwh.removeVertex)
 power                 rival-formula + planted-miss controls on tables recomputed at n = 13, 15, 17 (the ten fits minus the n = 19/21 ones, taken from E-299/E-301/E-304)
 orbits                n = 15 h = 3: are the image keys at K = 3 a function of the move orbit of the source row?
"""
import sys, time, collections, json
import numpy as np

def rows(m):
    L = m - 2; out = []
    def rec(i, last, cur):
        if i == L: out.append(tuple(cur)); return
        cur.append(0); rec(i+1, last, cur); cur.pop()
        for x in range(2, m - (i+1) + 1):
            e = i + 1 + x
            if e > last:
                cur.append(x); rec(i+1, e, cur); cur.pop()
    rec(0, 0, []); return out

def cartan(n, r):
    rel = [(i+1, i+1+x) for i, x in enumerate(r) if x]
    C = np.zeros((n, n))
    for s in range(1, n+1):
        stop = min((e for a, e in rel if a >= s), default=n+1)
        for t in range(s, min(n, stop-1)+1): C[t-1, s-1] = 1
    return C

def key(n, r):
    C = cartan(n, r); Phi = -np.linalg.solve(C.T, C); p = np.poly(Phi); q = np.rint(p)
    assert np.max(np.abs(p-q)) < 1e-6
    return tuple(int(v) for v in q)

def lone(n, h): return tuple([0]*h + [3] + [0]*(n-3-h))

def ktable(n, K1, hs, keep=False):
    """returns {h: (tab{(K,imgkeyindex):count}, rowsdata[(row,K,imgkey)] if keep)}"""
    tg = {}
    for h in hs: tg.setdefault(key(n, lone(n, h)), []).append(h)
    R = rows(n-K1); res = {k: [collections.Counter(), {}, []] for k in tg}
    for s in R:
        m = (0,)*K1 + s; k0 = key(n, m)
        if k0 not in tg: continue
        K = next((i for i, x in enumerate(m) if x), n)
        img = m[1:]; k = key(n-1, img); d = res[k0]; d[1].setdefault(k, len(d[1])); d[0][(K, d[1][k])] += 1
        if keep: d[2].append((m, K, d[1][k]))
    return {tuple(tg[k]): v for k, v in res.items()}

def splits(tab):
    """{K: number of image keys at that K}"""
    d = collections.defaultdict(set)
    for (K, j) in tab: d[K].add(j)
    return {K: len(v) for K, v in sorted(d.items())}

def K0of(tab, K1):
    """last K >= K1 with >= 2 image keys, or None"""
    s = splits(tab); c = [K for K, v in s.items() if v >= 2 and K >= K1]
    return max(c) if c else None

def show(n, K1, out):
    for hs, (tab, keys, _) in out.items():
        print(' n=%d class lone h in %s  rows %d  image keys %d' % (n, hs, sum(tab.values()), len(keys)))
        for K in sorted({k for k, _ in tab}):
            print('   K', K, {j: c for (kk, j), c in sorted(tab.items()) if kk == K})
        print('   K0 (last K>=%d with 2+ image keys) = %s ; law min(h,n-4-h) = %d' % (K1, K0of(tab, K1), min(min(hs), n-4-min(hs))))

mode = sys.argv[1]
if mode == 'table':
    n, K1 = int(sys.argv[2]), int(sys.argv[3]); hs = [int(x) for x in sys.argv[4].split(',')]; t = time.time()
    show(n, K1, ktable(n, K1, hs)); print('%.0fs' % (time.time()-t))
elif mode == 'power':
    # fits: (n, h, K1) ; h = short side. recorded n=19/21 fits are added as (n,h,K0) constants
    fits = {}
    for n, hs in ((13, (2, 3, 4)), (15, (3, 4, 5)), (17, (5, 6))):
        for h in hs:
            K1 = h - 1
            out = ktable(n, K1, [h])
            (tab, keys, _), = out.values()
            fits[(n, h)] = (K0of(tab, K1), tab, K1)
            print('refit n=%d h=%d K1=%d: K0=%s' % (n, h, K1, fits[(n, h)][0]), flush=True)
    rec = {(19, 5): 5, (19, 6): 6, (19, 7): 7, (21, 8): 8}   # recorded (E-299/E-301/E-304), central classes h=7 (n=19), h=8 (n=21)... see text
    law = lambda n, h: min(h, n-4-h)
    rivals = {'law min(h,n-4-h)': law,
              'R1 floor((n-1)/2)-2 (central formula for every class)': lambda n, h: (n-1)//2 - 2,
              'R2 min(h,n-5-h)': lambda n, h: min(h, n-5-h),
              'R3 min(h,n-3-h), small side': lambda n, h: min(h, n-3-h),
              'R3m min(h,n-3-h), read on long side h\'=n-4-h': lambda n, h: min(n-4-h, n-3-(n-4-h)),
              'R4 floor(n/2)-2': lambda n, h: n//2 - 2,
              'R5 h-1 (law minus one)': lambda n, h: h-1}
    obs = {k: v[0] for k, v in fits.items()}; obs.update(rec)
    even = {(10,2):2,(10,1):1,(12,3):3,(12,2):2,(14,4):4,(14,3):3,(14,2):2,(14,1):1,(16,5):5,(16,4):4}   # this round, experimentalist_even_extra_out.txt / n14_out.txt
    central_even = {(10,3):0,(12,4):0,(14,5):0}   # K0 = 0: no split at K >= 1
    obs_odd = dict(obs); obs.update(even)
    obs_c = dict(obs); obs_c.update(central_even)
    print('observed K0:', obs)
    for name, f in rivals.items():
        bad = [k for k, v in obs.items() if f(*k) != v]
        print('%-62s refuted on %d of %d odd-n fits %s' % (name, len(bad), len(obs_odd), [k for k in obs_odd if f(*k) != obs_odd[k]]))
        print('%-62s   +even non-central: refuted on %s of %d ; +central even: %s' % ('', [k for k in obs if f(*k) != obs[k]], len(obs), [k for k in central_even if f(*k) != central_even[k]]))
    # planted miss: instrument that misses the K0 split (image keys at K = K0 merged) must report K0-1 or less => law test says NO
    print('planted K0-1 miss (image keys of K = K0 merged):')
    for (n, h), (k0, tab, K1) in fits.items():
        t2 = collections.Counter()
        for (K, j), c in tab.items(): t2[(K, 0 if K == k0 else j)] += c
        k2 = K0of(t2, K1); print('  n=%d h=%d true K0 %s, planted-instrument K0 %s, law test %s' % (n, h, k0, k2, 'YES(pass)' if k2 == law(n, h) else 'NO(detected)'))
elif mode == 'orbits':
    from quivermutation import freeMoves, lnaMoves
    n, K1, h = 15, 2, 3; t = time.time()
    out = ktable(n, K1, [h], keep=True); (tab, keys, data), = out.values()
    S = {m: (K, j) for m, K, j in data}; print('class rows', len(S), '%.0fs' % (time.time()-t))
    parent = {m: m for m in S}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    outside = 0; edges = 0
    for m in S:
        for r in lnaMoves.rewritesOf(n, list(m)):
            r = tuple(r)
            if r in S:
                edges += 1; a, b = find(m), find(r)
                if a != b: parent[a] = b
            else: outside += 1
    print('rewrite edges inside class set', edges, 'outside (K<2 or not in class)', outside, '%.0fs' % (time.time()-t))
    comp = collections.defaultdict(list)
    for m in S: comp[find(m)].append(m)
    print('components of the K>=2 class rows under rewrites:', len(comp), 'sizes', sorted((len(v) for v in comp.values()), reverse=True)[:15])
    for K in (2, 3):
        cnt = collections.Counter()
        for r, v in comp.items():
            js = frozenset(S[m][1] for m in v if S[m][0] == K)
            if js: cnt[js] += 1
        print(' K=%d: components containing rows of this K, by set of image keys hit:' % K, dict(cnt))
        mixed = [v for v in comp.values() if len({S[m][1] for m in v if S[m][0] == K}) > 1]
        print('   components mixing both image keys at K=%d: %d' % (K, len(mixed)))
        if mixed:
            v = mixed[0]; print('   example component size', len(v), [m for m in v if S[m][0] == K][:3])
