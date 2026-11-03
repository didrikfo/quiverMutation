import sys,re,collections
sys.path.insert(0,"workshop/rounds/097")
from maverick_cartanZ import ct
import sympy as sp
from fractions import Fraction
def qhist(S):
    n=S.shape[0]; Si=S.inv(); cols=[[Fraction(int(Si[i,j].p),int(Si[i,j].q)) for i in range(n)] for j in range(n)]
    key=lambda v:tuple(x%1 for x in v)
    seen={key([Fraction(0)]*n)}; fr=list(seen)
    while fr:
        nf=[]
        for g in fr:
            for c in cols:
                h=key([a+b for a,b in zip(g,c)])
                if h not in seen: seen.add(h); nf.append(h)
        fr=nf
    Sl=[[int(S[i,j]) for j in range(n)] for i in range(n)]
    h=collections.Counter()
    for y in seen:
        q=sum(y[i]*Sl[i][j]*y[j] for i in range(n) for j in range(n))
        h[(q%2)]+=1   # y^T S y mod 2
    return len(seen),tuple(sorted(h.items()))
for f in ["skeptic_qsno.out","skeptic_qsno12.out"]:
  seen=set()
  for l in open("workshop/rounds/110/"+f):
    N=int(re.match(r"N=(\d+)",l).group(1)); a=eval(re.search(r"a=(\(.*?\)) b",l).group(1)); b=eval(re.search(r"b=(\(.*?\)) hist",l).group(1))
    if (N,a) in seen: continue
    seen.add((N,a))
    M=[sp.Matrix(__import__('numpy').array(ct.lnaCartanMatrix(N,list(r))).astype(int).tolist()) for r in (a,b)]
    S=[m+m.T for m in M]; ha,hb=qhist(S[0]),qhist(S[1])
    print(N,len(a),ha[0],hb[0],ha[1]==hb[1], ha[1],hb[1])
