"""R-terminal shapes: for short words, which placements are in the 444 orbit S (listed by right gap g = n - end of last relation, and left gap L = first start vertex - 1).
 Usage: timeout 10m .venv/bin/python workshop/rounds/019/theorist_shapes.py N"""
import sys
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
for w in ["4","45","46","47","48","44","334","335","357","3335","3345","3346","3347","4556","4667"]:
    os_=offs(w)
    ins=[(o,n-(o+len(w)-1+int(w[-1])+1)) for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
    print("n",n,w,"placements",len(os_),"in-S (offset L, right gap g):",ins)
