"""Double-mutation neighbours (interior, n = 20, offset 6) of a family of words, printed as word@offset.
  .venv/bin/python workshop/rounds/006/theorist_nbrs.py 34 4-8
"""
import sys
import batch
from quivermutation import doubleMutation
pre = sys.argv[1]; lo, hi = map(int, sys.argv[2].split("-")); n = 20; o = 6
def show(r):
    nz = [i for i, v in enumerate(r) if v]; a, b = nz[0], nz[-1]
    return "%s@%d" % ("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), a)
for x in range(lo, hi + 1):
    r = batch._rowFor(n, "%s%d" % (pre, x), o)
    print(pre + str(x), "@", o, "->", [show(v) + str(s[0]) for v, s in doubleMutation.rewritesOf(n, tuple(r))])
