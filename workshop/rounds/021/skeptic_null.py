"""Selectivity null for E-093/E-098. For n and each word length k in 4,5,6: ALL nondecreasing words over letters 2..9
(the only words that are LNAs), with >=4 placements. Per word: placements m, in-S count (S = 444 orbit rows), right gaps of in-S placements.
Strata: contains a 4 or not. Also binomial expectation of 'exactly one in S' from the stratum's per-placement in-S rate.
Usage: timeout 10m .venv/bin/python workshop/rounds/021/skeptic_null.py N   -> workshop/rounds/021/skeptic_null_n{N}.txt"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from collections import Counter, defaultdict
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
out=open("workshop/rounds/021/skeptic_null_n%d.txt"%n,"w")
def P(*a): print(*a,file=out,flush=True); print(*a,flush=True)
P("n",n,"|S|",len(S))
for k in (4,5,6):
    st=defaultdict(lambda: dict(words=0,pl=0,ins=0,one=0,oneg01=0,all=0,none=0,gaps=Counter(),exp=0.0,ex=[]))
    for t in itertools.combinations_with_replacement(range(2,10),k):
        w="".join(map(str,t)); os_=offs(w)
        if len(os_)<4: continue
        ins=[o for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
        key=("has4" if "4" in w else "no4")
        d=st[key]; d["words"]+=1; d["pl"]+=len(os_); d["ins"]+=len(ins)
        m=len(os_)
        if len(ins)==0: d["none"]+=1
        elif len(ins)==m: d["all"]+=1
        if len(ins)==1 and m>1:
            d["one"]+=1; g=n-(ins[0]+k+int(w[-1])); d["gaps"][g]+=1
            if g in (0,1): d["oneg01"]+=1
            d["ex"].append((w,ins[0],g))
        d.setdefault("ms",[]).append(m)
    for key,d in sorted(st.items()):
        p=d["ins"]/d["pl"]
        e=sum(m*p*(1-p)**(m-1) for m in d["ms"])
        P("k",k,key,"words",d["words"],"placements",d["pl"],"inS_rate %.3f"%p,"none",d["none"],"all",d["all"],
          "exactly_one",d["one"],"(binomial exp %.1f)"%e,"one&g<=1",d["oneg01"],"gaps",dict(sorted(d["gaps"].items())))
        P("   examples",d["ex"][:12])
