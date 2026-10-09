"""Which offsets of each PARTIAL word lie in S (the 444 orbit), n given. Usage: .venv/bin/python workshop/rounds/018/skeptic_partial_where.py N"""
import sys; sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves
n=int(sys.argv[1]); R=freeMoves.REDUCED
offs=lambda w:[o for o in range(n) if batch._rowFor(n,w,o) is not None]
o=offs("444"); S=frozenset(freeMoves.orbitReport(n,batch._rowFor(n,"444",o[len(o)//2]),free=R,limit=300000).rows)
for l in open("workshop/rounds/018/skeptic_rowset16_n%d_fast.txt"%n):
    f=l.split()
    if len(f)>3 and f[2]=="PARTIAL":
        w=f[0]; os_=offs(w)
        print(w,"offsets",os_[0],"..",os_[-1],"in S at",[x for x in os_ if freeMoves._startOf(batch._rowFor(n,w,x),R) in S])
