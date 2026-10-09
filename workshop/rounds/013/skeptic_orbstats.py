"""Stratified merge counts by word and by orbit from skeptic_orbscan_n{12..15}.txt (run from repo root)."""
from collections import defaultdict
def stratum(w):
    if "34" in w: return "has34"
    if "4" in w: return "4,no34"
    return "no4"
tot=defaultdict(lambda:[0,0,set(),set()]) # words, merged words, orbits(all), orbits(merged)
for n in (12,13,14,15):
    per=defaultdict(lambda:[0,0,set(),set()]); big=defaultdict(list)
    for l in open("workshop/rounds/013/skeptic_orbscan_n%d.txt"%n):
        w,no,m,h,oid,sz=l.split(); sz=int(sz)
        if sz==1: continue
        s=stratum(w); mg=(m=="merged")
        for d,key in ((per,(n,oid)),(tot,(n,oid))):
            d[s][0]+=1; d[s][1]+=mg; d[s][2].add(key)
            if mg: d[s][3].add(key)
        if mg: big[oid].append(w)
    print("n=%d"%n)
    for s in ("has34","4,no34","no4"):
        a=per[s]; print("  %-7s words %2d merged %2d | orbits hosting %2d, orbits hosting a merged word %2d"%(s,a[0],a[1],len(a[2]),len(a[3])))
    for o,ws in sorted(big.items(),key=lambda t:-len(t[1])): print("   merged orbit",o,len(ws),"words:",*ws)
print("pooled")
for s in ("has34","4,no34","no4"):
    a=tot[s]; print("  %-7s words %3d merged %3d | orbits %3d with merged %3d"%(s,a[0],a[1],len(a[2]),len(a[3])))
# extra: split no4 by presence of 2; and where the merged 4-words sit (orbit of 444 vs other)
print("extra")
ex=defaultdict(lambda:[0,0,set(),set()]); inbig=defaultdict(lambda:[0,0])
for n in (12,13,14,15):
    rows=[l.split() for l in open("workshop/rounds/013/skeptic_orbscan_n%d.txt"%n)]
    o444=[r[4] for r in rows if r[0]=="444"][0]
    for w,no,m,h,oid,sz in rows:
        if int(sz)==1: continue
        mg=m=="merged"
        s="no4,no2" if ("4" not in w and "2" not in w) else ("no4,has2" if "4" not in w else None)
        if s: ex[s][0]+=1; ex[s][1]+=mg; ex[s][2].add((n,oid)); ex[s][3].add((n,oid)) if mg else None
        if "4" in w and mg:
            k="has34" if "34" in w else "4,no34"
            inbig[(n,k)][0]+=1; inbig[(n,k)][1]+=(oid==o444)
for s,a in ex.items(): print("  %-9s words %3d merged %3d | orbits %3d with merged %3d"%(s,a[0],a[1],len(a[2]),len(a[3])))
for k,v in sorted(inbig.items()): print("  merged 4-words",k,"total",v[0],"in 444 orbit",v[1])
