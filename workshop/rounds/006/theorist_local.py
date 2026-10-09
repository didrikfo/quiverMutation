"""Is 33x@o ~ dual(33x)@(o+x-3) already reached by the FLOATING rules alone (no anchored rules, edge moves, double mutation)?
  timeout 10m .venv/bin/python workshop/rounds/006/theorist_local.py 20 5 5
args: n x offset.  Prints orbit size and whether the dual placement is in it, for floating-only and for the whole table.
"""
import sys
import batch
from quivermutation import freeMoves, lnaMoves
n, x, o = map(int, sys.argv[1:4])
R = freeMoves.REDUCED
a = batch._rowFor(n, "33%d" % x, o)
dual = freeMoves.mirrorRow(n, batch._rowFor(n, "33%d" % x, n - (x + 3) - (o + x - 3)))  # = dual core at offset o+x-3
want = freeMoves._startOf(dual, R)
for label, kw in (("floating+edges only", dict(rules=lnaMoves.VERIFIED_MOVES, edges=True, doubles=False)), ("floating+doubles only", dict(rules=lnaMoves.VERIFIED_MOVES, edges=False, doubles=True)), ("floating only", dict(rules=lnaMoves.VERIFIED_MOVES, edges=False, doubles=False)),
                  ("floating + edges + doubles", dict(rules=lnaMoves.VERIFIED_MOVES)),
                  ("whole table", dict())):
    w = freeMoves.orbitReport(n, a, free=R, limit=200000, **kw)
    print(label, "size", len(w.rows), w.stoppedBy, "dual@o+d in orbit:", want in w.rows,
          "c@s-o in orbit:", freeMoves._startOf(batch._rowFor(n, "33%d" % x, n - 2 * x - o), R) in w.rows)
