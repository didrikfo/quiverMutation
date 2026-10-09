"""For each split four-letter word (letters<=9, a 4, >=4 placements, exactly some placements in the 444 orbit S) at n:
find its in-S placement(s), BFS (reduced move set) to 333@0, print the path as words, and mark each step that is an R-step
(a,b,b,d)->(a-1,b,d+1) / (a,b,b)->(a-1,b) on the stripped word.  Also: positions of the in-S placement.
 Usage: timeout 10m .venv/bin/python workshop/rounds/019/theorist_split.py N"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from collections import deque
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
tgt=fm._startOf(tuple(batch._rowFor(n,"333",0)),R)
assert tgt in S
# BFS from target (moves are symmetric as a relation on S for connectivity; we use forward BFS from each source for honesty)
def show(r):
    nz=[i for i,x in enumerate(r) if x]
    if not nz: return "-"
    return "%s@%d"%("".join(str(x) if x<10 else "(%d)"%x for x in r[nz[0]:nz[-1]+1]),nz[0])
def bfs(src):
    par={src:None}; q=deque([src])
    while q:
        c=q.popleft()
        if c==tgt: break
        for m in fm.movesFrom(n,c,None,R):
            m=tuple(m)
            if m not in par: par[m]=c; q.append(m)
    if tgt not in par: return None
    p=[];c=tgt
    while c is not None: p.append(c); c=par[c]
    return p[::-1]
def isR(a,b):
    # a,b rows; is b the R-image of a at some position (compare as dense words)
    wa=[x for x in a]; 
    for i in range(len(wa)-2):
        x=wa[i:i+3]
        if x[1]==x[2] and x[1]>=3 and 3<=x[0]<=x[1] if False else False: pass
    return None
words=["".join(map(str,t)) for t in itertools.combinations_with_replacement(range(1,10),4) if 4 in t]
for w in words:
    os_=offs(w)
    if len(os_)<4: continue
    ins=[o for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
    if not ins or len(ins)==len(os_): continue
    print("WORD",w,"placements",len(os_),"0..%d"%(os_[-1]),"in S at",ins)
    for o in ins:
        p=bfs(fm._startOf(tuple(batch._rowFor(n,w,o)),R))
        print("  path len",None if p is None else len(p)-1,":"," -> ".join(show(r) for r in p) if p else "")
