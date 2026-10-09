import sys, time
sys.path.insert(0,'workshop'); 
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); R=freeMoves.REDUCED
for w in sys.argv[2:]:
    offs=[o for o in range(n) if batch._rowFor(n,w,o) is not None]
    t=time.time(); o=offs[len(offs)//2]
    rep=freeMoves.orbitReport(n,batch._rowFor(n,w,o),free=R,limit=300000)
    held=[p for p in offs if freeMoves._startOf(batch._rowFor(n,w,p),R) in rep.rows]
    print(n,w,"offs",len(offs),"from",o,"held",held,"size",len(rep.rows),rep.closed,round(time.time()-t,1),flush=True)
