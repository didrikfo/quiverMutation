"""Prediction test: the orbit of the single relation 3@2 (even n) / 3@3 (odd n) holds the word 35 (36) at exactly the
even offsets, and the orbit of 5@0 (6@0) holds it at exactly the odd offsets; the two are disjoint.
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_predict.py N [limit]
"""
import sys, time
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm
n = int(sys.argv[1]); limit = int(sys.argv[2]) if len(sys.argv) > 2 else 3000000
w, a, b = ('35', '3', '5') if n % 2 == 0 else ('36', '3', '6')
oa = 2 if n % 2 == 0 else 3
t = time.time()
res = {}
for name, (k, o) in (('P', (a, oa)), ('Q', (b, 0))):
    row = tuple(batch._rowFor(n, k, o))
    rep = fm.orbitReport(n, row, free=fm.REDUCED, limit=limit)
    held = [oo for oo in range(n) if batch._rowFor(n, w, oo) and fm._startOf(tuple(batch._rowFor(n, w, oo)), fm.REDUCED) in rep.rows]
    mir = fm._startOf(fm.mirrorRow(n, row), fm.REDUCED) in rep.rows
    res[name] = rep.rows
    print('n=%d orbit %s of %s@%d: size %d closed %s holds own mirror %s; %s at offsets %s (%.0fs)' % (n, name, k, o, len(rep.rows), rep.closed, mir, w, held, time.time() - t), flush=True)
print('disjoint:', not (res['P'] & res['Q']))
