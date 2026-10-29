"""Round 097: in the n = 15 key class of lone 3 at (5,6), which free ends with K >= 5 give which n = 14 image key (needs the hits files of experimentalist_scan15.py). Repo root."""
import glob, collections
from quivermutation import coxeterTables as ct, piecewiseHereditary as pwh
n = 15; D = 'workshop/rounds/097/'
rows = []
for f in sorted(glob.glob(D+'experimentalist_scan15_hits_5_*.txt')): rows += [tuple(int(c) for c in l.strip()) for l in open(f) if l.strip()]
keys = {}; tab = collections.Counter(); ex = {}
for m in rows:
    R = [(a+1, a+1+x) for a, x in enumerate(m) if x]
    Kh, Kt = min(s for s, e in R)-1, n-max(e for s, e in R)
    for side, K, v in (('head', Kh, 1), ('tail', Kt, n)):
        if K >= 5:
            k = ct.lnaCoxeterKey(n-1, tuple(pwh.removeVertex(n, list(m), v)[1])); j = keys.setdefault(k, len(keys))
            tab[(side, K, j, Kh >= 5 and Kt >= 5)] += 1; ex.setdefault((side, K, j), ''.join(map(str, m)))
for k, c in sorted(tab.items()): print(k, c)
print('examples', ex)
