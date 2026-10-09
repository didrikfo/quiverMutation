"""For 34x, 44x, 45x (x digit) at every offset: the orbit's held offsets of the same word, and the set of c = y + p of the 33y@p placements it contains
(33y conserves c = y + o, pairing c <-> n - c).  If c-sets are {c0, n-c0} the word is tied to the 33x chain by a fixed shift.
  timeout 10m .venv/bin/python workshop/rounds/006/theorist_link.py 14 34 4-9
"""
import sys
import batch
from quivermutation import freeMoves
n, pre = int(sys.argv[1]), sys.argv[2]
lo, hi = map(int, sys.argv[3].split("-"))
R = freeMoves.REDUCED
def cset(rows):
    cs = set()
    for r in rows:
        nz = [i for i, v in enumerate(r) if v]
        if len(nz) == 3 and nz[-1] - nz[0] == 2 and r[nz[0]] == 3 and r[nz[0]+1] == 3 and r[nz[2]] >= 3:
            cs.add(r[nz[2]] + nz[0])
    return sorted(cs)
for x in range(lo, hi + 1):
    word = "%s%d" % (pre, x)
    offs = [o for o in range(n) if batch._rowFor(n, word, o) is not None]
    done = set(); 
    for o in offs:
        if o in done: continue
        w = freeMoves.orbitReport(n, batch._rowFor(n, word, o), free=R, limit=300000)
        held = [p for p in offs if freeMoves._startOf(batch._rowFor(n, word, p), R) in w.rows]
        done.update(held)
        print(n, word, "hi", offs[-1], "offset", o, "held", held, "c-set of 33y", cset(w.rows), "closed" if w.closed else "CAP", len(w.rows), flush=True)
