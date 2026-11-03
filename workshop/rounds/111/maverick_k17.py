"""S-1 / H-026, round 111: K-threshold law at n = 17, key level, WITHOUT the full class scan.
Idea: a row with a free end of K >= K1 vertices is 0^K1 + (row of length n-K1) (head) or (row of length n-K1) + 0^K1 (tail);
so all rows with an end K >= K1 are enumerated from allRelationLengths(n-K1). Filter by key == key(lone 3 at (h,K)), h+K=n-4.
usage (repo root): maverick_k17.py check            # shift hypothesis vs full enumeration at n=12 (K1=3,4)
                   maverick_k17.py run N K1         # K-table of head/tail-end deletion images, key class of the lone 3 with h=(N-4)//2 (lower h)
Key level only (lnaCoxeterKey); images by pwh.removeVertex at v=1 / v=N."""
import sys, collections, time
from quivermutation import coxeterTables as ct, nakayama, piecewiseHereditary as pwh
def lone(n, h): return tuple([0]*h + [3] + [0]*(n-2-h-1))
mode = sys.argv[1]
if mode == 'check':
    n = 12
    full = [tuple(r) for r in nakayama.allRelationLengths(n)]
    for K1 in (3, 4):
        sub = [tuple(r) for r in nakayama.allRelationLengths(n-K1)]
        head = {tuple(r) for r in full if all(x == 0 for x in r[:K1])}
        tail = {tuple(r) for r in full if all(x == 0 for x in r[len(r)-K1:])}
        print('n', n, 'K1', K1, 'full head', len(head), 'shift', len(sub), 'equal', head == {(0,)*K1+s for s in sub},
              '| tail', len(tail), 'equal', tail == {s+(0,)*K1 for s in sub})
else:
    n, K1 = int(sys.argv[2]), int(sys.argv[3]); t0 = time.time()
    h0 = (n-4)//2
    tg = {ct.lnaCoxeterKey(n, lone(n, h)) for h in (h0, n-4-h0)}
    print('target keys', len(tg), 'lone rows', [(h, n-4-h) for h in range(n-3) if ct.lnaCoxeterKey(n, lone(n, h)) in tg])
    sub = [tuple(r) for r in nakayama.allRelationLengths(n-K1)]
    rows = {(0,)*K1+s for s in sub}   # head-end rows only; tail ends are their mirrors (class mirror-closed; tail trailing-zero rows are NOT a shift, see check)
    print('rows with an end >= %d: %d (%.0fs)' % (K1, len(rows), time.time()-t0))
    cls = [r for r in rows if ct.lnaCoxeterKey(n, r) in tg]
    print('in key class:', len(cls), '(%.0fs)' % (time.time()-t0))
    tab = collections.Counter(); keys = {}; none = 0
    for m in cls:
        R = [(a+1, a+1+x) for a, x in enumerate(m) if x]
        ends = (min(s for s, e in R)-1, n-max(e for s, e in R)) if R else None
        for K, v in zip(ends[:1], (1,)):
            if K >= K1:
                x = pwh.removeVertex(n, list(m), v)
                if x is None: none += 1; continue
                k = ct.lnaCoxeterKey(n-1, tuple(x[1])); keys.setdefault(k, len(keys)); tab[(K, keys[k])] += 1
    for K in range(K1, n):
        d = {j: c for (kk, j), c in sorted(tab.items()) if kk == K}
        if d: print(' K', K, d)
    print(' None:', none)
    for K0 in range(K1, n): print(' image keys for K>=%d: %d' % (K0, len({j for (kk, j) in tab if kk >= K0})))
    print('%.0fs' % (time.time()-t0))
