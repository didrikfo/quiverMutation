"""S-1 item 3, round 097: n = 15 key scan for the lone-3 key classes (5,6)/(6,5) and (4,7)/(7,4), K-threshold K0 = 5 (E-146, E-173).
usage (repo root): experimentalist_scan15.py --plan | scan SHARD NSHARDS | ends
scan: rows with index i % NSHARDS == SHARD; key of every row (no deletions); hits -> workshop/rounds/097/experimentalist_scan15_hits_<h>_<shard>.txt
ends: union of hits, class sizes, free-end (K >= 3..) images by removeVertex at n = 14, key table K x image-key. Key level only."""
import sys, glob, collections, time
from quivermutation import coxeterTables as ct, nakayama, piecewiseHereditary as pwh
n = 15
D = 'workshop/rounds/097/'
def lone(h): return tuple([0]*h + [3] + [0]*(n-2-h-1))
HS = (4, 5, 6, 7)   # (h,K) with h+K = n-4 = 11: (4,7) (5,6) (6,5) (7,4)
TG = {h: ct.lnaCoxeterKey(n, lone(h)) for h in HS}
mode = sys.argv[1]
if mode == '--plan':
    print('target keys equal: (5,6)=(6,5)', TG[5] == TG[6], '(4,7)=(7,4)', TG[4] == TG[7], '(4,7)=(5,6)', TG[4] == TG[5])
    t = time.time(); c = 0
    for lna in nakayama.allRelationLengths(n):
        c += 1
        if c == 3000: break
    print('3000 rows enumerated in %.2fs' % (time.time()-t))
    t = time.time(); c = 0
    for lna in nakayama.allRelationLengths(n):
        c += 1
        ct.lnaCoxeterKey(n, tuple(lna))
        if c == 3000: break
    print('3000 keys (first rows) %.2fs' % (time.time()-t))
elif mode == 'scan':
    sh, ns = int(sys.argv[2]), int(sys.argv[3]); t = time.time(); hits = collections.defaultdict(list); c = 0
    kk = collections.defaultdict(list)
    for h, v in TG.items(): kk[v].append(h)
    for lna in nakayama.allRelationLengths(n):
        c += 1
        if c % ns != sh: continue
        k = ct.lnaCoxeterKey(n, tuple(lna))
        if k in kk: hits[kk[k][0]].append(lna)
    print('shard', sh, 'scanned', c, {h: len(v) for h, v in hits.items()}, '%.0fs' % (time.time()-t))
    for h in HS: open(D+'experimentalist_scan15_hits_%d_%d.txt' % (h, sh), 'w').write('\n'.join(''.join(map(str, r)) for r in hits.get(h, [])))
else:
    for h in (5, 4):
        rows = []
        for f in sorted(glob.glob(D+'experimentalist_scan15_hits_%d_*.txt' % h)):
            rows += [tuple(int(c) for c in l.strip()) for l in open(f) if l.strip()]
        print('=== key class of lone 3 at h =', h, 'size', len(rows), 'lone in class', [x for x in HS if lone(x) in set(rows)])
        tab = collections.Counter(); keys = {}; nonecount = collections.Counter()
        for m in rows:
            R = [(a+1, a+1+x) for a, x in enumerate(m) if x]
            ends = (min(s for s, e in R)-1, n-max(e for s, e in R)) if R else None
            for K, v in zip(ends, (1, n)):
                if K >= 3:
                    x = pwh.removeVertex(n, list(m), v)
                    if x is None: nonecount[(h, min(K, 7))] += 1; continue
                    im = tuple(x[1]); k = ct.lnaCoxeterKey(n-1, im)
                    keys.setdefault(k, len(keys)); tab[(min(K, 7), keys[k])] += 1
        for K in (3, 4, 5, 6, 7): print(' K', K if K < 7 else '>=7', {j: c for (kk_, j), c in sorted(tab.items()) if kk_ == K})
        print(' removeVertex returned None (skipped):', dict(nonecount))
        for K0 in (3, 4, 5, 6): print(' image keys for K>=%d:' % K0, len({j for (kk_, j) in tab if kk_ >= K0}))
