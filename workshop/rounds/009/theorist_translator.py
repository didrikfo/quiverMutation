"""Closure criterion for drift families: an orbit is 'offset-rigid' (its own word is held at fewer than all offsets)
iff it holds no TRANSLATOR, a short word (span <= 3 nonzero-support letters) held at EVERY offset at which the word has a row.
For each prefix word p+x (x digit) and each offset o, compute the orbit (reduced walk), list translators held, own-word offsets held.
  timeout 10m .venv/bin/python workshop/rounds/009/theorist_translator.py 14 35 5-9
"""
import sys
import batch
from quivermutation import freeMoves
n, pre = int(sys.argv[1]), sys.argv[2]
lo, hi = map(int, sys.argv[3].split("-"))
R = freeMoves.REDUCED
def wordsOf(rows, maxspan=3):
    out = {}
    for r in rows:
        nz = [i for i, v in enumerate(r) if v]
        if not nz: continue
        a, b = nz[0], nz[-1]
        if b - a + 1 <= maxspan:
            out.setdefault("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), set()).add(a)
    return out
_offs = {}
def offsOf(word):
    if word not in _offs:
        _offs[word] = {o for o in range(n) if batch._rowFor(n, word, o) is not None}
    return _offs[word]
def translators(rows):
    return sorted(w for w, s in wordsOf(rows).items() if "(" not in w and len(offsOf(w)) >= 4 and s >= offsOf(w))
agree = tot = 0
for x in range(lo, hi + 1):
    word = "%s%d" % (pre, x)
    offs = sorted(offsOf(word)); done = set()
    for o in offs:
        if o in done: continue
        w = freeMoves.orbitReport(n, batch._rowFor(n, word, o), free=R, limit=300000)
        held = [p for p in offs if freeMoves._startOf(batch._rowFor(n, word, p), R) in w.rows]
        done.update(held)
        T = translators(w.rows)
        merged = set(held) == set(offs)
        print(n, word, "o", o, "held", held, "merged_all" if merged else "rigid", "T", T, "closed" if w.closed else "CAP", len(w.rows), flush=True)
        if len(offs) >= 4:
            tot += 1; agree += (bool(T) == merged)
print("agree translator<=>all-offsets-merged", agree, "/", tot)
