"""Orbit scan of the 4-letter words of the --max-word 4 catalogue that contain a 4 (letters 1..9, nondecreasing, >= MINOFF placements) at length n.
Usage (repo root): .venv/bin/python workshop/rounds/014/experimentalist_orbscan4.py N [MINOFF] [LIMIT]
Output workshop/rounds/014/experimentalist_orbscan4_n{N}.txt, line: word noffs merged|rigid nheld orbit_id orbit_size closed
"""
import sys, itertools
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); MINOFF=int(sys.argv[2]) if len(sys.argv)>2 else 4
LIMIT=int(sys.argv[3]) if len(sys.argv)>3 else 300000
R=freeMoves.REDUCED
out=open("workshop/rounds/014/experimentalist_orbscan4_n%d.txt"%n,"w")
sizes=[]; seen={}
def oid_of(row):
    st=freeMoves._startOf(row,R)
    if st not in seen:
        rep=freeMoves.orbitReport(n,row,free=R,limit=LIMIT)
        k=len(sizes); sizes.append((len(rep.rows),rep.closed))
        for s in rep.rows: seen[s]=k
        seen[st]=k
    return seen[st]
for t in itertools.combinations_with_replacement(range(1,10),4):
    if 4 not in t: continue
    w="".join(map(str,t))
    offs=[o for o in range(n) if batch._rowFor(n,w,o) is not None]
    if len(offs)<MINOFF: continue
    ids=[oid_of(batch._rowFor(n,w,o)) for o in offs]
    mid=ids[len(offs)//2]
    held=sum(1 for i in ids if i==mid)
    print(w,len(offs),"merged" if held==len(offs) else "rigid",held,mid,sizes[mid][0],int(sizes[mid][1]),file=out,flush=True)
print("done",file=out)
