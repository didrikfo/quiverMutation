"""4046-type cores (4046 5046 5056) at one n: held classes by orbit, sizes; tests a period-2 translation (offset classes by parity, steps of 2).
  timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_core4046.py N WORD [LIMIT]   (WORD like 4046 or 504)"""
import sys, time
import batch
from quivermutation import freeMoves
n, word = int(sys.argv[1]), sys.argv[2]
limit = int(sys.argv[3]) if len(sys.argv) > 3 else 300000
R = freeMoves.REDUCED
offs = [o for o in range(n) if batch._rowFor(n, word, o) is not None]
done, out = set(), []
for o in offs:
    if o in done: continue
    t0 = time.time()
    w = freeMoves.orbitReport(n, batch._rowFor(n, word, o), free=R, limit=limit)
    held = [p for p in offs if freeMoves._startOf(batch._rowFor(n, word, p), R) in w.rows]
    done.update(held)
    print(n, word, "offset", o, "held", held, "size", len(w.rows), "closed" if w.closed else "CAP", "%.0fs" % (time.time()-t0), flush=True)
