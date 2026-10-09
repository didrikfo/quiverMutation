"""Small orbit (235/255/455/2455/3334) vs 444 orbit: list the short words (<=4 relations, gaps allowed) held in each, and a labelled shortest path 3334 -> 235 and 2455 -> 235.
  timeout 10m .venv/bin/python workshop/rounds/017/theorist_small.py N
"""
import sys
sys.path.insert(0,'.'); sys.path.insert(0,'workshop/rounds/013')
import batch
from collections import Counter
from quivermutation import freeMoves as fm
from theorist_path import bfs
n=int(sys.argv[1]); R=fm.REDUCED
def orb(w):
    for o in range(n):
        r=batch._rowFor(n,w,o)
        if r is not None: break
    return fm.orbitReport(n,tuple(r),free=R,limit=300000)
def word(r): return ''.join(str(x) for x in r if x)
def ct(r): return sum(1 for x in r if x)
for name,w in (("small","235"),("big","444")):
    rep=orb(w); print(name,"size",len(rep.rows),rep.closed)
    c=Counter(ct(r) for r in rep.rows); print(" #relations histogram",sorted(c.items()))
    ws=Counter(word(r) for r in rep.rows if ct(r)<=4)
    print(" words with <=4 relations (word:count):",' '.join("%s:%d"%kv for kv in sorted(ws.items(),key=lambda kv:(len(kv[0]),kv[0]))))
    # which rows have every entry gap-free (contiguous word)?
for a,b in (("3334","235"),("2455","235"),("3334","2455")):
    oa=[o for o in range(n) if batch._rowFor(n,a,o) is not None]; ob=[o for o in range(n) if batch._rowFor(n,b,o) is not None]
    res=bfs(n,tuple(batch._rowFor(n,a,oa[len(oa)//2])),tuple(batch._rowFor(n,b,ob[len(ob)//2])))
    print("PATH",a,oa[len(oa)//2],"->",b,ob[len(ob)//2])
    if isinstance(res,str) or res is None: print(res); continue
    s,p=res; print(' ',''.join(map(str,s)))
    for r,l in p: print(' ',''.join(map(str,r)),'<-',l)
