"""Orbit partition, key partition and orbit sizes of one word's offsets at length n (full walk of each orbit).
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_word.py N word [word..]"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm, coxeterTables as ct
n = int(sys.argv[1])
for w in sys.argv[2:]:
    orb = {}; keys = {}; sizes = []
    for o in range(n):
        row = batch._rowFor(n, w, o)
        if row is None: continue
        s = fm._startOf(tuple(row), fm.REDUCED)
        if s not in orb:
            rep = fm.orbitReport(n, tuple(row), free=fm.REDUCED, limit=1500000)
            i = len(sizes); sizes.append((len(rep.rows), rep.closed))
            for r in rep.rows: orb[r] = i
        keys.setdefault(ct.lnaCoxeterKey(n, tuple(row)), []).append(o)
    byorb = {}
    for o in range(n):
        row = batch._rowFor(n, w, o)
        if row is not None: byorb.setdefault(orb[fm._startOf(tuple(row), fm.REDUCED)], []).append(o)
    print('n=%d %s orbits(offsets,size,closed): %s | key classes: %s' % (n, w, [(v, sizes[k]) for k, v in byorb.items()], list(keys.values())), flush=True)
