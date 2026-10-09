"""One-step neighbours of 33x@o under the double mutation alone, and under the floating table alone (no free-move additions).
  .venv/bin/python workshop/rounds/006/theorist_step.py 18 5 6
"""
import sys
import batch
from quivermutation import doubleMutation, lnaMoves
n, x, o = map(int, sys.argv[1:4])
def show(r):
    r = list(r); a = next((i for i, v in enumerate(r) if v), None)
    b = max(i for i, v in enumerate(r) if v)
    return "%s@%d" % ("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), a)
row = batch._rowFor(n, "33%d" % x, o)
print("row", show(row))
for r, seq in doubleMutation.rewritesOf(n, tuple(row)): print("  double", show(r), seq)
for r in lnaMoves.rewritesOf(n, list(row), lnaMoves.VERIFIED_MOVES): print("  table ", show(r))
