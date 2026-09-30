"""Single-step neighbours (double mutation) of 2-letter words at interior offset: which xy slide (xy@o -> xy@o+1 in one or two steps)?
  .venv/bin/python workshop/rounds/009/theorist_translate.py
"""
import batch
from quivermutation import doubleMutation
n, o = 20, 6
def show(r):
    nz = [i for i, v in enumerate(r) if v]
    if not nz: return "0"
    a, b = nz[0], nz[-1]
    return "%s@%d" % ("".join(str(v) if v < 10 else "(%d)" % v for v in r[a:b+1]), a)
for w in ["23", "34", "44", "45", "56", "24", "35"]:
    r = batch._rowFor(n, w, o)
    if r is None: print(w, "no row"); continue
    print(w, sorted(show(v) for v, s in (doubleMutation.rewritesOf(n, tuple(r)) or [])))
