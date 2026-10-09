"""Like round 010 skeptic_scan.py but records the orbit id of each word (orbit at the middle offset) and the orbit id of every held offset.
Usage (repo root): python workshop/rounds/013/skeptic_orbscan.py N [MINOFF]  -> workshop/rounds/013/skeptic_orbscan_n{N}.txt
Line: word noffs merged|rigid nheld orbit_id orbit_size
"""
import sys, itertools, time
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); MINOFF=int(sys.argv[2]) if len(sys.argv)>2 else 4
R=freeMoves.REDUCED
out=open("workshop/rounds/013/skeptic_orbscan_n%d.txt"%n,"w")
sizes=[]; seen={}
def oid_of(row):
    st=freeMoves._startOf(row,R)
    if st not in seen:
        rep=freeMoves.orbitReport(n,row,free=R,limit=300000)
        assert rep.closed
        k=len(sizes); sizes.append(len(rep.rows))
        for s in rep.rows: seen[s]=k
        seen[st]=k
    return seen[st]
for x,y,z in itertools.product(range(1,10),repeat=3):
    w="%d%d%d"%(x,y,z)
    offs=[o for o in range(n) if batch._rowFor(n,w,o) is not None]
    if len(offs)<MINOFF: continue
    ids=[oid_of(batch._rowFor(n,w,o)) for o in offs]
    mid=ids[len(offs)//2]
    held=sum(1 for i in ids if i==mid)
    print(w,len(offs),"merged" if held==len(offs) else "rigid",held,mid,sizes[mid],file=out,flush=True)
