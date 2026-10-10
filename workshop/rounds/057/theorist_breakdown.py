"""T2 (referee item 5): trimmed-word histogram of the rows of the 444 orbit, n given. Row = tuple of relation lengths;
word = row with leading/trailing zeros removed (rows with a part >= 10 shown as lists). Also: 33y rows with y+p in {3,n-3}?"""
import sys, collections
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); R=freeMoves.REDUCED
o4=[o for o in range(n) if batch._rowFor(n,"444",o) is not None]
r=freeMoves.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000); assert r.closed
H=collections.Counter(); ex={}
for row in r.rows:
    t=list(row)
    while t and t[0]==0: t.pop(0)
    while t and t[-1]==0: t.pop()
    w="".join(map(str,t)) if all(v<10 for v in t) else str(t)
    H[w]+=1; ex.setdefault(w,row)
print("n",n,"rows",len(r.rows),"distinct words",len(H))
for w,c in sorted(H.items(),key=lambda x:(len(x[0]),x[0])): 
    if w.startswith("33") and len(w)==3 or len(w)<=3 or c>0 and len(w)<=4: print(w,c)
print("words of length>=5:",sum(c for w,c in H.items() if len(w)>=5))
