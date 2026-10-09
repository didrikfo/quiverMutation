"""For each split 4-letter word (letters<=9, a 4, >=4 placements, partly in the 444 orbit S) at n: take its in-S placement(s) and ask
whether repeated application of lemma R alone (plus deleting 2-arrow relations) reaches 333@0.  R: intervals at consecutive starts
s,s+1,s+2[,s+3] of lengths a,b,b[,d], 3<=a<=b, -> (s,s+a-1),(s+1,s+1+b)[,(s+2,s+3+d)].  Also R-chain for ALL placements of the word.
 Usage: timeout 10m .venv/bin/python workshop/rounds/019/theorist_rchain.py N"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from collections import deque
from quivermutation import freeMoves as fm, doubleMutation as dm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
tgt=fm._startOf(tuple(batch._rowFor(n,"333",0)),R)
def show(r):
    nz=[i for i,x in enumerate(r) if x]
    return "-" if not nz else "%s@%d"%("".join(str(x) if x<10 else "(%d)"%x for x in r[nz[0]:nz[-1]+1]),nz[0])
def rmoves(r):
    iv=dm.intervalsOf(r); st={a:(a,b) for a,b in iv}; out=[]
    for (s,e) in iv:
        a=e-s; i1=st.get(s+1); i2=st.get(s+2)
        if not(i1 and i2): continue
        b=i1[1]-i1[0]
        if i2[1]-i2[0]!=b or not 3<=a<=b: continue
        i3=st.get(s+3)
        new=set(iv)-{(s,e),i1,i2}|{(s,s+a-1),(s+1,s+1+b)}
        if i3:
            new-= {i3}; new|={(s+2,s+3+(i3[1]-i3[0]))}
        res=dm.relLengthsOf(n,new)
        if res is not None and fm.stripLengthTwo(res) in {fm.stripLengthTwo(tuple(v)) for v,_ in dm.rewritesOf(n,r)}: out.append(fm.stripLengthTwo(res))
    return out
def chain(src):
    par={src:None}; q=deque([src])
    while q:
        c=q.popleft()
        if c==tgt: break
        for m in rmoves(c):
            if m not in par: par[m]=c; q.append(m)
    if tgt in par:
        p=[];c=tgt
        while c is not None: p.append(c); c=par[c]
        return p[::-1],None
    term=[r for r in par if not rmoves(r)]
    return None,term,len(par)
words=["".join(map(str,t)) for t in itertools.combinations_with_replacement(range(1,10),4) if 4 in t]
print("n",n,"|S|",len(S))
tot=ok=0
for w in words:
    os_=offs(w)
    if len(os_)<4: continue
    ins=[o for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
    if not ins or len(ins)==len(os_): continue
    for o in ins:
        src=fm._startOf(tuple(batch._rowFor(n,w,o)),R); tot+=1
        p,*rest=chain(src)
        if p: ok+=1
        print(w,"placement %d of 0..%d"%(o,os_[-1]),"R-chain to 333@0:",(" -> ".join(show(r) for r in p)) if p else "NO; R-closure size %d, R-terminals %s (in S: %s)"%(rest[1],[show(t) for t in rest[0]],[t in S for t in rest[0]]))
print("split words with in-S placement reducing by R to 333@0: %d of %d"%(ok,tot))
