"""Which core placements (word, offset) lie in the reduced-walk orbit of word@offset?  Only words of <= maxlen letters, printed grouped by word with (offset) lists.
  timeout 10m .venv/bin/python workshop/rounds/006/theorist_words.py 14 345 3 [maxlen]
"""
import sys
import batch
from quivermutation import freeMoves
n, word, o = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]); maxlen = int(sys.argv[4]) if len(sys.argv) > 4 else 4
R = freeMoves.REDUCED
w = freeMoves.orbitReport(n, batch._rowFor(n, word, o), free=R, limit=300000)
out = {}
for r in w.rows:
    nz = [i for i, v in enumerate(r) if v]
    if not nz: continue
    a, b = nz[0], nz[-1]
    if b - a + 1 <= maxlen:
        out.setdefault("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), []).append(a)
print(n, word, o, "orbit", len(w.rows), "closed" if w.closed else "CAP")
for k in sorted(out, key=lambda k: (len(k), k)): print(" ", k, sorted(out[k]))
