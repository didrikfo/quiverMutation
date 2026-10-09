"""Position null: for the exactly-one-in-S words of skeptic_null.py (re-derived), the expected chance that a UNIFORMLY chosen placement has right gap <=1
 (sum over words of #placements with g<=1 / m) vs observed count. Also exact binomial/ Poisson-binomial tail P(X>=obs).
Usage: timeout 10m .venv/bin/python workshop/rounds/021/skeptic_null_gap.py N"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
from collections import defaultdict
agg=defaultdict(lambda:[0,0,0.0,[]])
for k in (4,5,6):
    for t in itertools.combinations_with_replacement(range(2,10),k):
        w="".join(map(str,t)); os_=offs(w)
        if len(os_)<4: continue
        ins=[o for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
        if len(ins)!=1: continue
        g=lambda o:n-(o+k+int(w[-1]))
        key=(k,"has4" if "4" in w else "no4"); a=agg[key]
        p=sum(1 for o in os_ if g(o)<=1)/len(os_)
        a[0]+=1; a[1]+= g(ins[0])<=1; a[2]+=p; a[3].append(p)
def tail(ps,x):
    dp=[1.0]
    for p in ps:
        nd=[0.0]*(len(dp)+1)
        for i,v in enumerate(dp): nd[i]+=v*(1-p); nd[i+1]+=v*p
        dp=nd
    return sum(dp[x:])
for key,a in sorted(agg.items()):
    print("n",n,key,"exactly-one words",a[0],"with g<=1:",a[1],"expected under uniform placement %.2f"%a[2],"P(X>=obs) %.3f"%tail(a[3],a[1]))
