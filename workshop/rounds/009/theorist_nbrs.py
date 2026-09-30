"""One-step neighbours of word@o (interior, n=20, o=6) under (a) the double mutation, (b) the reduced walk (table+free+edges) minus (a).
  .venv/bin/python workshop/rounds/009/theorist_nbrs.py 44 4-7
"""
import sys
import batch
from quivermutation import doubleMutation, freeMoves
pre = sys.argv[1]; lo, hi = map(int, sys.argv[2].split("-")); n = 20; o = 6
def show(r):
    r = list(r); nz = [i for i, v in enumerate(r) if v]
    if not nz: return "0"
    a, b = nz[0], nz[-1]
    return "%s@%d" % ("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), a)
for x in range(lo, hi + 1):
    r = batch._rowFor(n, "%s%d" % (pre, x), o)
    dm = {show(v) for v, s in doubleMutation.rewritesOf(n, tuple(r))}
    rw = freeMoves.reducedMovesFrom(n, tuple(r))
    rw = {show(v[0] if isinstance(v, tuple) and v and isinstance(v[0], (tuple, list)) else v) for v in rw}
    print(pre + str(x), "@", o, "double:", sorted(dm), "| reduced-walk other:", sorted(rw - dm))
