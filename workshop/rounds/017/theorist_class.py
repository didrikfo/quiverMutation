"""Class label j of the orbit of each placement of each word: j = offsets o' of 333@o' in the orbit (set; {} = none). Odd n (no parity split).
  timeout 10m .venv/bin/python workshop/rounds/017/theorist_class.py N WORD...
"""
import sys
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
t333={o:fm._startOf(tuple(batch._rowFor(n,"333",o)),R) for o in offs("333")}
cache={}
def orb(row):
    st=fm._startOf(tuple(row),R)
    if st in cache: return cache[st]
    rep=fm.orbitReport(n,tuple(row),free=R,limit=400000); assert rep.closed
    v=(len(rep.rows),tuple(sorted(o for o,s in t333.items() if s in rep.rows)))
    for s in rep.rows: cache[s]=v
    return v
for w in sys.argv[2:]:
    ds={}
    for o in offs(w):
        sz,lab=orb(batch._rowFor(n,w,o)); ds.setdefault((sz,lab),[]).append(o)
    print(w," | ".join("size %d class%s at %s"%(sz,list(lab),v) for (sz,lab),v in ds.items()))
