"""Does k@0 lie in the reduced orbit of 3@(k-2)? (relation of k arrows at offset 0 vs. relation of 3 arrows ending at the same
vertex minus...)  Prints for n, k: same orbit?  Also 3@(k-3).   .venv/bin/python workshop/rounds/011/theorist_k0.py nmin nmax
"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    out = []
    for k in range(4, 10):
        r0 = batch._rowFor(n, str(k), 0)
        if r0 is None: continue
        rep = fm.orbitReport(n, tuple(r0), free=fm.REDUCED, limit=300000)
        if not rep.closed: out.append('%d:capped' % k); continue
        hit = [j for j in range(0, n) if batch._rowFor(n, '3', j) and fm._startOf(tuple(batch._rowFor(n, '3', j)), fm.REDUCED) in rep.rows]
        out.append('%d@0 holds 3@%s (size %d)' % (k, hit, len(rep.rows)))
    print('n=%d: %s' % (n, '; '.join(out)), flush=True)
