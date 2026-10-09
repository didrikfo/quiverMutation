"""Shortest reduced-walk path between 33x at offset o and at offset s-o (s = n-2x), with the rule used per step.
  .venv/bin/python workshop/rounds/006/theorist_path.py 13 4 0      # n, x, offset
"""
import sys
from collections import deque
import batch
from quivermutation import freeMoves, lnaMoves
n, x, o = map(int, sys.argv[1:4])
word = "33%d" % x
s = n - 2 * x
R = freeMoves.REDUCED
def row(off): return freeMoves._startOf(batch._rowFor(n, word, off), R)
a, b = row(o), row(s - o)
print("start", o, a); print("target", s - o, b)
par = {a: None}; q = deque([a])
while q and b not in par:
    r = q.popleft()
    for m in freeMoves.reducedMovesFrom(n, r):
        if m not in par:
            par[m] = r; q.append(m)
    if len(par) > 400000: break
if b not in par: print("not found", len(par)); sys.exit()
path = []; r = b
while r is not None: path.append(r); r = par[r]
path.reverse()
print("len", len(path) - 1, "orbit seen", len(par))
for r in path: print("".join(str(c) if c < 10 else "(%d)" % c for c in r))
