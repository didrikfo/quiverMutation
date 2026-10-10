"""T1: is k(c) (E-062 / round-004 shift table) a simple function of the core word? Leave-one-out over 11 cores
(k from experimentalist_shift_table.txt n=13, plus 33x family k=2x from E-065)."""
import itertools, re, numpy as np
rows={}
for l in open('workshop/rounds/004/experimentalist_shift_table.txt'):
    p=l.split()
    if p[1]=='13' and p[2]=='s': rows[p[0]]=int(p[7])
for x in range(3,9): rows.setdefault('33%d'%x,2*x)
print(len(rows),rows)
def feats(w):
    d=[int(c) for c in w]; nz=[v for v in d if v]
    return dict(one=1,len=len(d),sum=sum(d),first=d[0],last=d[-1],mx=max(d),mn=min(nz),nz=len(nz),zeros=len(d)-len(nz),
                lastm1=d[-1]-d[-2] if len(d)>1 else 0)
names=list(feats('33').keys()); names.remove('one')
import sys
if len(sys.argv)>1 and sys.argv[1]=="no33": rows={w:k for w,k in rows.items() if not w.startswith("33") or len(w)==4}; rows.pop("334",None)
ws=sorted(rows); y=np.array([rows[w] for w in ws],float)
X={w:feats(w) for w in ws}
res=[]
for r in (1,2,3):
    for S in itertools.combinations(names,r):
        A=np.array([[1]+[X[w][f] for f in S] for w in ws],float)
        err=0;ex=0
        for i in range(len(ws)):
            m=np.arange(len(ws))!=i
            c=np.linalg.lstsq(A[m],y[m],rcond=None)[0]
            e=A[i]@c-y[i]; err+=e*e; ex+=abs(round(A[i]@c)-y[i])<1e-9
        res.append((err,ex,S))
res.sort()
for r in res[:6]: print(r)
print('best exact LOO hits out of',len(ws),':',max(r[1] for r in res))
