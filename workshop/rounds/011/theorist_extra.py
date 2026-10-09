"""Orbits (size, own offsets) and keys of given words at n, to see the parity-class orbits not in P, Q.
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_extra.py N word [word..]"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm, coxeterTables as ct
n = int(sys.argv[1])
for w in sys.argv[2:]:
    done = {}
    for o in range(n):
        row = batch._rowFor(n, w, o)
        if row is None: continue
        s = fm._startOf(tuple(row), fm.REDUCED)
        if s in done: print(w, o, 'orbit', done[s]); continue
        rep = fm.orbitReport(n, tuple(row), free=fm.REDUCED, limit=3000000)
        for r in rep.rows: done[r] = (len(done) and max(v[0] for v in done.values()) + 1 or 0, len(rep.rows))
        done[s] = (max(v[0] for v in done.values()), len(rep.rows))
        sing = [(k, oo) for k in range(3, n - 1) for oo in range(n) if batch._rowFor(n, str(k), oo) and fm._startOf(tuple(batch._rowFor(n, str(k), oo)), fm.REDUCED) in rep.rows]
        print(w, o, 'orbit', done[s], 'closed', rep.closed, 'single relations in it:', sing)
