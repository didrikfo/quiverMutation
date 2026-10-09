"""Orbit id (by size, with 444 orbit marked) of every single-relation row a@o, and of the words 35, 455, 3334, 2455, 444, 333 placements, n given.
  timeout 10m .venv/bin/python workshop/rounds/017/theorist_single.py N
"""
import sys
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
cache={}
def oid(row):
    st=fm._startOf(tuple(row),R)
    if st in cache: return cache[st]
    rep=fm.orbitReport(n,tuple(row),free=R,limit=400000)
    assert rep.closed
    i=(len(rep.rows),min(rep.rows))
    for s in rep.rows: cache[s]=i
    return i
big=oid(batch._rowFor(n,"444",[o for o in range(n) if batch._rowFor(n,"444",o) is not None][0]))
ids={big:"BIG"}
def lab(i):
    if i not in ids: ids[i]="O%d(%d)"%(len(ids),i[0])
    return ids[i]
lab(big)
for w in ["3","4","5","6","7","8","9","10","35","455","3334","2455","333","444","36","34","45","55","33","44"]:
    out=[]
    for o in range(n):
        r=batch._rowFor(n,w,o)
        if r is not None: out.append("%d:%s"%(o,lab(oid(r))))
    print(w,' '.join(out))
print(ids)
