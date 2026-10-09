"""Closure test for a drift family P x.  Seed = P x0 (x0 = last letter of P).  Collapse products = span<=2 double-mutation/reduced neighbours of the seed.
PREDICTION: family P x (x >= x0) is offset-merged iff for SOME INTERIOR offset o the orbit of seed@o holds a collapse product of seed@o at every offset where that word has a row
(a 'full translator reached by collapse').  OUTCOME: the orbit of P x@o holds P x@o' for all offsets o' (merged) or not (rigid, chain-paired).
  timeout 10m .venv/bin/python workshop/rounds/009/theorist_closure.py 14 55 5-8
"""
import sys
import batch
from quivermutation import doubleMutation, freeMoves
n, pre = int(sys.argv[1]), sys.argv[2]
lo, hi = map(int, sys.argv[3].split("-"))
R = freeMoves.REDUCED
def show(r):
    nz = [i for i, v in enumerate(r) if v]
    if not nz: return None, None
    a, b = nz[0], nz[-1]
    return "".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), a
def offsOf(word):
    return {o for o in range(n) if batch._rowFor(n, word, o) is not None}
def wordsHeld(rows, span):
    out = {}
    for r in rows:
        w, a = show(r)
        if w and len(w) <= span: out.setdefault(w, set()).add(a)
    return out
def orbit(word, o):
    return freeMoves.orbitReport(n, batch._rowFor(n, word, o), free=R, limit=300000)
seed = pre + pre[-1]
print("n", n, "seed", seed)
pred = False
for o in sorted(offsOf(seed)):
    r = batch._rowFor(n, seed, o)
    prods = set()
    for v, s in (doubleMutation.rewritesOf(n, tuple(r)) or []):
        w, a = show(v)
        if w and len(w) <= 2 and "(" not in w: prods.add(w)
    w_ = orbit(seed, o); held = wordsHeld(w_.rows, 2)
    full = sorted(p for p in prods if held.get(p, set()) >= offsOf(p) and len(offsOf(p)) >= 3)
    print(" seed@%d" % o, "products", sorted(prods), "full translators among them", full, "orbit", len(w_.rows), "closed" if w_.closed else "CAP", flush=True)
    if 0 < o < max(offsOf(seed)): pred = pred or bool(full)   # interior offsets only: edge-touching seeds have extra moves
print("PREDICT merged" if pred else "PREDICT rigid")
for x in range(lo, hi + 1):
    word = "%s%d" % (pre, x); offs = sorted(offsOf(word)); done = set(); res = []
    for o in offs:
        if o in done: continue
        w = orbit(word, o)
        held = [p for p in offs if freeMoves._startOf(batch._rowFor(n, word, p), R) in w.rows]
        done.update(held); res.append((held, len(w.rows), w.closed))
    print(word, "OUTCOME", "merged" if len(res) == 1 and len(offs) > 1 else "rigid", [(h, s, "closed" if c else "CAP") for h, s, c in res], flush=True)
