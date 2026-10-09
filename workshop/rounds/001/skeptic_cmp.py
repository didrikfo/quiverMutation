import sys,time
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm, coxeterTables as ct
def part(n,core,limit=1500000):
    offs=[o for o in range(n) if batch._rowFor(n,core,o)]
    st={o:fm._startOf(tuple(batch._rowFor(n,core,o)),fm.REDUCED) for o in offs}
    reps={}
    for o in offs:
        reps[o]=fm.orbitReport(n,tuple(batch._rowFor(n,core,o)),free=fm.REDUCED,limit=limit)
    orb=[];seen=set()
    for o in offs:
        if o in seen: continue
        g=[p for p in offs if st[p] in reps[o].rows]; seen|=set(g); orb.append((g,len(reps[o].rows),reps[o].stoppedBy))
    keys={}
    for o in offs: keys.setdefault(ct.lnaCoxeterKey(n,tuple(batch._rowFor(n,core,o))),[]).append(o)
    return orb,list(keys.values())
n=int(sys.argv[1])
for core in sys.argv[2:]:
    t=time.time();o,k=part(n,core)
    print(n,core,'orbits',o,'coxeter',k,'same' if sorted(g for g,_,_ in o)==sorted(k) else 'DIFFER',round(time.time()-t))
    sys.stdout.flush()
