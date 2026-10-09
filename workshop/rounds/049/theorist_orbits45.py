"""T4: orbit sizes of core `45` at every offset, reduced walk, to see which offsets share an orbit (F-051/F-053 mechanism).
usage: theorist_orbits45.py N [limit]"""
import sys
from quivermutation import lnaMoves as lm, freeMoves as fm
n = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
for o in range(0, n - 5):
    row = lm.embedPattern(n, ((0, 4), (1, 5)), o) if False else None
    rl = [0]*(n-2); 
    if o + 1 >= len(rl): break
    rl[o] = 4; rl[o+1] = 5
    if not lm.isAdmissible(n, rl): print(o, 'not admissible'); continue
    w = fm.orbitReport(n, rl, free=fm.REDUCED, limit=lim)
    print(o, w.stoppedBy, len(w.rows), hash(frozenset(w.rows)) % 100000, flush=True)
