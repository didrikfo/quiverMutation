"""Merged-ness of every 3-letter word xyz (letters 1..9) at one n: 'merged' = the reduced-walk orbit of the word at one interior offset holds the word at ALL its offsets (>= MINOFF offsets).
Orbits are cached by row set so each orbit is walked once.  Usage: skeptic_scan.py N [MINOFF] [letters-lo-hi]   -> skeptic_scan_n{N}.txt
"""
import sys, itertools, time
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); MINOFF=int(sys.argv[2]) if len(sys.argv)>2 else 4
lo,hi=map(int,(sys.argv[3] if len(sys.argv)>3 else "1-9").split("-"))
R=freeMoves.REDUCED
out=open("workshop/rounds/010/skeptic_scan_n%d.txt"%n,"w")
done_orbits=[]
seen={}   # start -> (id,size,closed)
t0=time.time()
for x,y,z in itertools.product(range(lo,hi+1),repeat=3):
    w="%d%d%d"%(x,y,z)
    offs=[o for o in range(n) if batch._rowFor(n,w,o) is not None]
    if len(offs)<MINOFF: continue
    o=offs[len(offs)//2]; st=freeMoves._startOf(batch._rowFor(n,w,o),R)
    if st not in seen:
        rep=freeMoves.orbitReport(n,batch._rowFor(n,w,o),free=R,limit=300000)
        k=len(done_orbits); done_orbits.append(1)
        for s in rep.rows: seen[s]=(k,len(rep.rows),rep.closed)
        seen[st]=(k,len(rep.rows),rep.closed)
    oid,size,closed=seen[st]
    held=[p for p in offs if freeMoves._startOf(batch._rowFor(n,w,p),R) in seen and seen[freeMoves._startOf(batch._rowFor(n,w,p),R)][0]==oid]
    print(w,len(offs),"merged" if len(held)==len(offs) else "rigid",len(held),size,"closed" if closed else "CAP",file=out,flush=True)
print("done",round(time.time()-t0),file=out,flush=True)
