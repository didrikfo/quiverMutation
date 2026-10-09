"""Size and time of the reduced orbit of one word placement (capped).  timeout 10m .venv/bin/python workshop/rounds/013/theorist_size.py N WORD OFFSET LIMIT"""
import sys, time
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm
n, w, o, lim = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
t = time.time()
rep = fm.orbitReport(n, tuple(batch._rowFor(n, w, o)), free=fm.REDUCED, limit=lim)
print('n=%d %s@%d: %d rows, %s, %.0fs' % (n, w, o, len(rep.rows), rep.stoppedBy, time.time() - t))
