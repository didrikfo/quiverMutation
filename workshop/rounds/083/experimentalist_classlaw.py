"""S-1 class-level K-threshold law (round 083). usage: experimentalist_classlaw.py N
Key class of the lone 3 at the middle (h,K) (h+K = N-4) of N-vertex LNAs, scanned over all rows by lnaCoxeterKey (as E-146).
Free ends (K = free vertices at head/tail end, as maverick_n12) with K >= 3 are deleted (removeVertex); the image is labelled by
its derivedOrbits(N-1) orbit (class-level, a lower bound for the derived class: different orbits may be one class) and by its key.
Prints, for each threshold k, the number of distinct image orbits / keys over ends with K >= k."""
import sys, time, collections
from quivermutation import coxeterTables as ct, nakayama, piecewiseHereditary as pwh, freeMoves
n = int(sys.argv[1]); t0 = time.time()
tg = {}
for h in range(0, n - 3):
    tg[ct.lnaCoxeterKey(n, tuple([0]*h + [3] + [0]*(n-h-1)))] = tg.get(ct.lnaCoxeterKey(n, tuple([0]*h + [3] + [0]*(n-h-1))), []) + [h]
rows = {k: [] for k in tg}
for lna in nakayama.allRelationLengths(n):
    k = ct.lnaCoxeterKey(n, tuple(lna))
    if k in rows: rows[k].append(tuple(lna))
print('n', n, 'scan %.0fs' % (time.time() - t0), 'key classes of lone 3:', {tuple(v): len(rows[k]) for k, v in tg.items()}, flush=True)
_, orbs = freeMoves.derivedOrbits(n - 1)
lab = {r: i for i, v in enumerate(orbs.values()) for r in v}
EXT = '--ext' in sys.argv
if EXT:
    _, orbs2 = freeMoves.derivedOrbits(n - 1, edges=True, doubles=True)
    lab2 = {r: i for i, v in enumerate(orbs2.values()) for r in v}
    print('ext orbits at n-1:', len(orbs2), 'plain', len(orbs), flush=True)
print('derivedOrbits(%d) %.0fs' % (n - 1, time.time() - t0), flush=True)
for k, hs in tg.items():
    rs = rows[k]
    if len(rs) < 2: continue
    # orbit structure at level n (class level) of this key class
    ext = collections.defaultdict(set); tab = collections.defaultdict(collections.Counter); kt = collections.defaultdict(set)
    for m in rs:
        R = [(a+1, a+1+x) for a, x in enumerate(m) if x]
        if not R: continue
        ends = (min(s for s, e in R) - 1, n - max(e for s, e in R))
        for K, v in zip(ends, (1, n)):
            if K >= 2:
                im = tuple(pwh.removeVertex(n, list(m), v)[1]); kk = min(K, 7)
                tab[kk][lab[im]] += 1; ext[kk].add(lab2[im]) if EXT else None; kt[kk].add(ct.lnaCoxeterKey(n - 1, im))
    print('key class lone-3 at h in', hs, 'rows', len(rs))
    for thr in range(2, 8):
        orb = set(); ks = set(); cnt = 0; eo = set()
        for kk in range(thr, 8): eo |= ext[kk]; orb |= set(tab[kk]); ks |= kt[kk]; cnt += sum(tab[kk].values())
        if cnt: print('  ends K >=', thr, ': ends', cnt, 'image orbits', len(orb), 'image keys', len(ks), ('ext orbits %d' % len(eo)) if EXT else '')
print('done %.0fs' % (time.time() - t0))
