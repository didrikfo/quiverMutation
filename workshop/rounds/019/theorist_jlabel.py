"""For each split 4-letter word at n: for every placement o, the class label j = min of offsets o' with 333@o' in its orbit (orbit size too).
Prints the in-S offset (j=0), the end gap g = n - (end vertex of the last relation) at the in-S placement, and whether j is a linear function of o.
 Usage: timeout 10m .venv/bin/python workshop/rounds/019/theorist_jlabel.py N"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
t333={o:fm._startOf(tuple(batch._rowFor(n,"333",o)),R) for o in offs("333")}
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
cache={}
def lab(row):
    st=fm._startOf(tuple(row),R)
    if st in cache: return cache[st]
    rep=fm.orbitReport(n,tuple(row),free=R,limit=60000)
    if not rep.closed: v=("cap",0)
    else:
        js=sorted(o for o,s in t333.items() if s in rep.rows)
        v=(min(js) if js else -1,len(rep.rows))
    if rep.closed:
        for s in rep.rows: cache[s]=v
    cache[st]=v
    return v
words=["".join(map(str,t)) for t in itertools.combinations_with_replacement(range(1,10),4) if 4 in t]
for w in words:
    os_=offs(w)
    if len(os_)<4: continue
    ins=[o for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
    if not ins or len(ins)==len(os_): continue
    omax=os_[-1]; o=ins[0]
    row=batch._rowFor(n,w,o); last=max(i+1+o+int(c) for i,c in enumerate(w))  # end vertex of the last relation (1-based starts)
    labs=[lab(batch._rowFor(n,w,x)) for x in os_]
    print(w,"omax",omax,"in-S o*",o,"o*-omax",o-omax,"gap",n-last,"| j@o:"," ".join("%d:%s/%s"%(x,l[0],l[1]) for x,l in zip(os_,labs)),flush=True)
