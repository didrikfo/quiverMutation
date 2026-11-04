"""Round 113, experimentalist: independent re-run of the S-1 / H-026 K-threshold (E-296) at key level, head ends, lone-3 class.
Different code path from rounds/111/maverick_k17.py: own row enumerator (DFS on the antichain condition, no `lines`/`nakayama`),
own Cartan builder, own Coxeter polynomial (numpy: charpoly of -C^{-T}C, rounded, residual checked), own vertex-1 deletion
(drop the first entry of the row). Cross-checks against the library are in `xcheck`.
usage (repo root):
  experimentalist_kthresh.py plan N K1        # row count of the shifted enumeration, no key computation
  experimentalist_kthresh.py xcheck           # own enumerator/key/deletion vs library at n=9..11
  experimentalist_kthresh.py run N K1 [libdel]  # K-table of image keys; libdel also compares with pwh.removeVertex on class rows
"""
import sys, time, collections
import numpy as np

def rows(m):
    """all relation-length rows of an LNA with m vertices: m-2 entries, x=0 or x>=2, start i (1-based) ends i+x <= m, ends strictly increasing."""
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
    C = cartan(n, r)
    Phi = -np.linalg.solve(C.T, C)
    p = np.poly(Phi)
    q = np.rint(p)
    assert np.max(np.abs(p-q)) < 1e-6, np.max(np.abs(p-q))
    return tuple(int(v) for v in q)

def lone(n, h): return tuple([0]*h + [3] + [0]*(n-3-h))

mode = sys.argv[1]
if mode == 'plan':
    n, K1 = int(sys.argv[2]), int(sys.argv[3]); t = time.time(); R = rows(n-K1)
    print('shifted rows at length', n-K1, ':', len(R), 'enum %.1fs' % (time.time()-t))
    t = time.time(); [key(n, (0,)*K1+r) for r in R[:300]]; print('300 keys %.2fs -> est total %.0fs' % (time.time()-t, (time.time()-t)*len(R)/300))
elif mode == 'xcheck':
    from quivermutation import nakayama, coxeterTables as ct, piecewiseHereditary as pwh
    import random
    for m in range(5, 12):
        a = set(rows(m)); b = {tuple(x) for x in nakayama.allRelationLengths(m)}
        print('enum m', m, len(a), len(b), a == b)
    random.seed(1); bad = 0; N = 0
    for n in (9, 10, 11):
        R = rows(n)
        for r in random.sample(R, min(400, len(R))):
            N += 1; bad += key(n, r) != ct.lnaCoxeterKey(n, r)
    print('key vs library on', N, 'rows, mismatches', bad)
    bad = N = 0
    for n in (9, 10):
        for r in rows(n):
            if r[0] == 0:
                x = pwh.removeVertex(n, list(r), 1); N += 1; bad += (x is None or tuple(x[1]) != r[1:])
    print('delete vertex 1: own (drop first entry) vs library on', N, 'head-free rows, mismatches', bad)
else:
    n, K1 = int(sys.argv[2]), int(sys.argv[3]); libdel = len(sys.argv) > 4; t0 = time.time()
    h0 = (n-4)//2
    tg = {key(n, lone(n, h)) for h in (h0, n-4-h0)}
    print('target keys', len(tg), 'lone rows', [(h, n-4-h) for h in range(n-3) if key(n, lone(n, h)) in tg])
    R = rows(n-K1); print('shifted rows', len(R), '%.0fs' % (time.time()-t0))
    tab = collections.Counter(); keys = {}; ncls = 0; mism = 0
    for s in R:
        m = (0,)*K1 + s
        if key(n, m) not in tg: continue
        ncls += 1
        K = next((i for i, x in enumerate(m) if x), None)
        img = m[1:]
        if libdel:
            from quivermutation import piecewiseHereditary as pwh
            x = pwh.removeVertex(n, list(m), 1); mism += (x is None or tuple(x[1]) != img)
        k = key(n-1, img); keys.setdefault(k, len(keys)); tab[(K, keys[k])] += 1
    print('in key class', ncls, '(%.0fs)' % (time.time()-t0))
    for K in range(K1, n):
        d = {j: c for (kk, j), c in sorted(tab.items()) if kk == K}
        if d: print(' K', K, d)
    for K0 in range(K1, n):
        v = {j for (kk, j) in tab if kk >= K0}
        if v: print(' image keys for K>=%d: %d' % (K0, len(v)))
    if libdel: print(' libdel mismatches', mism)
    print('%.0fs' % (time.time()-t0))
