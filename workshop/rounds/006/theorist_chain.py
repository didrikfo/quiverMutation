"""Check the chain picture for 33x: (a) double-mutation neighbours of 33x@o include 33(x-1)@(o+1) and 33(x+1)@(o-1) (interior);
(b) orbit(33x@o) holds the whole chain S_c = {33x'@(c-x')}, c = x+o, and 33x'@(n-c-x') for x' <= n-c, nothing else of the form 33y@p.
  timeout 10m .venv/bin/python workshop/rounds/006/theorist_chain.py 14
"""
import sys
import batch
from quivermutation import doubleMutation, freeMoves
n = int(sys.argv[1]); R = freeMoves.REDUCED
def red(r): return freeMoves._startOf(r, R)
def row(x, o): return batch._rowFor(n, "33%d" % x, o)
bad = 0; tot = 0
for x in range(4, min(n - 4, 10)):
    for o in range(0, n - x - 3 + 1):
        r = row(x, o)
        if r is None: continue
        nb = {tuple(v) for v, _ in doubleMutation.rewritesOf(n, tuple(r))}
        lo = row(x - 1, o + 1); hi = row(x + 1, o - 1) if o >= 1 else None
        okd = (tuple(lo) in nb) and (hi is None or tuple(hi) in nb or row(x + 1, o - 1) is None)
        if not okd: print("double-mutation step missing", x, o, (tuple(lo) in nb), hi is None or tuple(hi) in nb)
        tot += 1; bad += (not okd)
print("double-mutation steps checked", tot, "missing", bad)
# orbit membership
for x in range(3, 9):
    hi_off = n - x - 3
    if hi_off < 0: continue
    for o in range(0, hi_off + 1):
        w = freeMoves.orbitReport(n, row(x, o), free=R, limit=300000)
        c = x + o
        held = []  # every (x', o') with 33x'@o' in the orbit
        for xp in range(3, min(n - 3, 10)):
            for op in range(0, n - xp - 2):
                r = row(xp, op)
                if r is not None and red(r) in w.rows: held.append((xp, op))
        chain = [(xp, c - xp) for xp in range(3, min(c, 9) + 1) if c - xp <= n - xp - 3]
        pred = set(chain) | {(xp, n - c - xp) for xp in range(3, min(n - c, 9) + 1) if 0 <= n - c - xp <= n - xp - 3}
        print(n, x, o, "closed" if w.closed else "CAP", "size", len(w.rows), "held==pred" if set(held) == pred else "DIFF %s vs %s" % (sorted(held), sorted(pred)))
