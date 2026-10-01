"""Single-relation rows (k, o) whose Coxeter key equals the key of 35@0 (even n) or 36@0 (odd n), n = 8..20.
Keys only: no orbit walks.   .venv/bin/python workshop/rounds/011/theorist_triples.py [nmin nmax]
Also: the key of every placement of the word (35 or 36) and whether it is that key at every offset.
"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import coxeterTables as ct
lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (8, 20)
for n in range(lo, hi + 1):
    w = '35' if n % 2 == 0 else '36'
    ref = ct.lnaCoxeterKey(n, tuple(batch._rowFor(n, w, 0)))
    allsame = all(ct.lnaCoxeterKey(n, tuple(batch._rowFor(n, w, o))) == ref for o in range(n) if batch._rowFor(n, w, o))
    hits = []
    for k in range(2, n - 1):
        for o in range(n):
            r = batch._rowFor(n, str(k), o)
            if r is not None and ct.lnaCoxeterKey(n, tuple(r)) == ref: hits.append((k, o))
    ks = sorted(set(k for k, o in hits))
    print('n=%d word %s all offsets one key: %s; single relations (k,o): %s' % (n, w, allsame, hits), flush=True)
