"""(1) one-step double-mutation neighbours: is 35 (stripped) a neighbour of every placement of 3334 and 2455 (and which offset)?
(2) c-label of each placement's orbit: the offsets o' with 333@o' in the orbit (c = o'+3 convention below: c(333@o')=o'+3) and orbit size.
  timeout 10m .venv/bin/python workshop/rounds/017/theorist_label.py N
"""
import sys
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm, doubleMutation
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
print("n",n)
# (1)
for w,tgt in (("3334","35"),("2455","35"),("455","35")):
    res=[]
    for o in offs(w):
        r=tuple(batch._rowFor(n,w,o))
        nb={fm.stripLengthTwo(tuple(v)) for v,_ in doubleMutation.rewritesOf(n,r)}
        hit=[o2 for o2 in offs(tgt) if fm.stripLengthTwo(tuple(batch._rowFor(n,tgt,o2))) in nb]
        res.append((o,hit))
    print("double-mutation nbr",w,"->",tgt,res)
# (2)
t333={o:fm._startOf(tuple(batch._rowFor(n,"333",o)),R) for o in offs("333")}
cache={}
def orb(row):
    st=fm._startOf(tuple(row),R)
    if st in cache: return cache[st]
    rep=fm.orbitReport(n,tuple(row),free=R,limit=400000); assert rep.closed
    v=(len(rep.rows),tuple(sorted(o for o,s in t333.items() if s in rep.rows)))
    for s in rep.rows: cache[s]=v
    return v
for w in ("35","455","3334","2455","3335","444","34","334","36","333"):
    parts=[]
    for o in offs(w):
        sz,lab=orb(batch._rowFor(n,w,o))
        parts.append("%d:%d{%s}"%(o,sz,",".join(map(str,lab))))
    print(w," ".join(parts))
