"""Shortest reduced-walk path (n = 14) from 34@5 to 44@5, and the list of 2-letter words xy (2<=x,y<=9) whose double-mutation neighbours include xy shifted by 1.
  .venv/bin/python workshop/rounds/009/theorist_path34.py
"""
import batch
from collections import deque
from quivermutation import doubleMutation, freeMoves
R = freeMoves.REDUCED
def show(r):
    nz = [i for i, v in enumerate(r) if v]
    if not nz: return "0"
    a, b = nz[0], nz[-1]
    return "%s@%d" % ("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), a)
n = 20; o = 6
for x in range(2, 10):
    for y in range(2, 10):
        r = batch._rowFor(n, "%d%d" % (x, y), o)
        if r is None: continue
        nb = [show(v) for v, s in (doubleMutation.rewritesOf(n, tuple(r)) or [])]
        if ("%d%d@%d" % (x, y, o + 1)) in nb or ("%d%d@%d" % (x, y, o - 1)) in nb: print("self-sliding 2-letter word", x, y)
n = 14
src = tuple(batch._rowFor(n, "34", 5)); tgt = tuple(batch._rowFor(n, "44", 5))
src = tuple(freeMoves._startOf(src, R)); tgt = tuple(freeMoves._startOf(tgt, R))
par = {src: None}; q = deque([src])
while q:
    u = q.popleft()
    if u == tgt: break
    for v in freeMoves.reducedMovesFrom(n, u):
        v = tuple(v[0]) if v and isinstance(v[0], (tuple, list)) else tuple(v)
        if v not in par: par[v] = u; q.append(v)
path = []; u = tgt
while u is not None and u in par: path.append(show(u)); u = par[u]
print(" <- ".join(path))
