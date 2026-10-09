"""S-1 at n = 12: key classes of the lone 3 at (h,K)=(4,4) and (3,5); free-end (K>=3 / K>=4) images at n = 11 by key; move orbits in the class.
usage: python workshop/rounds/042/maverick_n12.py --plan | scan shard nshards | ends   (hits in maverick_n12class_<shard>.txt)"""
import sys, glob, collections, time
from quivermutation import coxeterTables as ct, freeMoves, lnaMoves, edgeMoves, doubleMutation, nakayama, piecewiseHereditary as pwh
n = 12
D = 'workshop/rounds/042/'
def lone(h): return tuple([0]*h + [3] + [0]*(n-2-h-1))
TG = {h: ct.lnaCoxeterKey(n, tuple([0]*h+[3]+[0]*(n-h-1))) for h in (3, 4, 5)}
mode = sys.argv[1]
if mode == '--plan':
    t = time.time(); c = 0
    for lna in nakayama.allRelationLengths(n):
        c += 1
        if c == 2000: break
    print('rows total 58786 (counted); 2000 keys in %.1fs' % (time.time()-t), 'target keys equal (4,4) vs (3,5):', TG[4] == TG[3], 'vs (5,3):', TG[3]==TG[5])
elif mode == 'scan':
    sh, ns = int(sys.argv[2]), int(sys.argv[3]); t = time.time(); hits = collections.defaultdict(list); c = 0
    kk = {v: h for h, v in TG.items()}
    for lna in nakayama.allRelationLengths(n):
        c += 1
        if c % ns != sh: continue
        k = ct.lnaCoxeterKey(n, tuple(lna))
        if k in kk: hits[kk[k]].append(lna)
    print('scanned', c, {h: len(v) for h, v in hits.items()}, '%.0fs' % (time.time()-t))
    for h, v in hits.items(): open(D+'maverick_n12class_%d_%d.txt' % (h, sh), 'w').write('\n'.join(''.join(map(str, r)) for r in v))
else:
    for h in (4, 5):
        rows = []
        for f in sorted(glob.glob(D+'maverick_n12class_%d_*.txt' % h)):
            rows += [tuple(int(c) for c in l.strip()) for l in open(f) if l.strip()]
        print('=== key class of lone 3 at h =', h, 'size', len(rows), 'lone in class', lone(h) in rows, 'lone(3,5)', lone(3) in rows, lone(5) in rows)
        idx = {m: i for i, m in enumerate(rows)}; par = list(range(len(rows)))
        def find(a):
            while par[a] != a: par[a] = par[par[a]]; a = par[a]
            return a
        for m in rows:
            cand = list(lnaMoves.rewritesOf(n, list(m), lnaMoves.ALL_MOVES)); cand.append(freeMoves.stripLengthTwo(m))
            cand += [r for r, _ in edgeMoves.rewritesOf(n, m)] + [r for r, _ in doubleMutation.rewritesOf(n, m)]
            for r in cand:
                r = tuple(r)
                if r in idx:
                    a, b = find(idx[m]), find(idx[r])
                    if a != b: par[a] = b
        orb = collections.Counter(find(i) for i in range(len(rows))); print('orbits', sorted(orb.values(), reverse=True))
        for h2 in ((4,) if h == 4 else (3, 5)):
            if lone(h2) in idx: print(' lone 3 h=%d orbit size %d' % (h2, orb[find(idx[lone(h2)])]))
        tab = collections.Counter(); keys = {}
        for i, m in enumerate(rows):
            R = [(a+1, a+1+x) for a, x in enumerate(m) if x]
            ends = (min(s for s, e in R)-1, n-max(e for s, e in R)) if R else None
            for K, v in zip(ends, (1, n)):
                if K >= 3:
                    im = tuple(pwh.removeVertex(n, list(m), v)[1]); k = ct.lnaCoxeterKey(n-1, im)
                    keys.setdefault(k, len(keys)); tab[(min(K, 6), keys[k])] += 1
        for K in (3, 4, 5, 6): print(' K', K if K < 6 else '>=6', {j: c for (kk, j), c in sorted(tab.items()) if kk == K})
        print(' image keys for K>=4:', sorted({j for (kk, j) in tab if kk >= 4}), ' for K>=3:', sorted({j for (kk, j) in tab}))
