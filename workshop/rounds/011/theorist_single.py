"""Orbit sizes and Coxeter keys of the single-relation rows `k@o` at length n (reduced walk).
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_single.py N
"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm, coxeterTables as ct
n = int(sys.argv[1]); seen = {}
keys = {}
for k in range(3, n - 1):
    for o in range(n):
        row = batch._rowFor(n, str(k), o)
        if row is None: continue
        s = fm._startOf(tuple(row), fm.REDUCED)
        key = ct.lnaCoxeterKey(n, tuple(row))
        kid = keys.setdefault(key, len(keys))
        if s not in seen:
            rep = fm.orbitReport(n, tuple(row), free=fm.REDUCED, limit=300000)
            oid = len(set(v[0] for v in seen.values()))
            for r in rep.rows: seen[r] = (oid, len(rep.rows), rep.closed)
        print('%s@%d orbit %d size %d closed %s key %d' % (k, o, seen[s][0], seen[s][1], seen[s][2], kid), flush=True)
