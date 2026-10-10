"""T2: is the 333@0 orbit (size 3767 at n=14) the same row set as the 444 orbit (E-083/E-091 size 3767)? usage: N"""
import sys
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); R=freeMoves.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
def orb(w,o):
    r=freeMoves.orbitReport(n,batch._rowFor(n,w,o),free=R,limit=300000); assert r.closed; return frozenset(r.rows)
o4=offs("444"); S=orb("444",o4[len(o4)//2])
print("n",n,"|444 orbit|",len(S))
for o in offs("333"):
    T=orb("333",o); print("333@%d"%o,len(T),"equal444" if T==S else ("subset" if T<S else "overlap %d"%len(T&S)))
