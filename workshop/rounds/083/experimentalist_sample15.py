"""Per-row cost of the derivedOrbits loop body (rewritesOf + stripLengthTwo) at n = 14 and 15 on S evenly spaced rows (default 3000), to size derivedOrbits(15).
usage: experimentalist_sample15.py [S]"""
import sys, time, math
from quivermutation import nakayama, lnaMoves, freeMoves
S = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
for n in (12, 13, 14, 15):
    R = math.comb(2*(n-1), n-1)//n; st = max(1, R//S); rows = []
    t = time.time()
    for i, r in enumerate(nakayama.allRelationLengths(n)):
        if i % st == 0: rows.append(r)
    te = time.time() - t; t = time.time(); nm = 0
    for r in rows:
        nm += len(list(lnaMoves.rewritesOf(n, list(r), lnaMoves.ALL_MOVES))); freeMoves.stripLengthTwo(r)
    tm = (time.time() - t) / len(rows)
    print('n', n, 'rows', R, 'enum %.1fs' % te, 'sampled', len(rows), 'moves/row %.1f' % (nm/len(rows)), 'loop body %.3f ms/row -> %.0f s' % (tm*1e3, tm*R), flush=True)
