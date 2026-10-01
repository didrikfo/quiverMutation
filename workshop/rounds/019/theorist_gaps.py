"""Membership only: for each split 4-letter word at n, the in-S placement o*, left gap L=o* (first relation starts at vertex L+1) and right gap g = n - end of last relation.
 Usage: timeout 10m .venv/bin/python workshop/rounds/019/theorist_gaps.py N   -> prints 'word L g' per split word, plus whether the word is in-S only at (L,g)"""
import sys, itertools
sys.path.insert(0,'.')
import batch
from quivermutation import freeMoves as fm
n=int(sys.argv[1]); R=fm.REDUCED
def offs(w): return [o for o in range(n) if batch._rowFor(n,w,o) is not None]
o4=offs("444"); S=frozenset(fm.orbitReport(n,batch._rowFor(n,"444",o4[len(o4)//2]),free=R,limit=300000).rows)
words=["".join(map(str,t)) for t in itertools.combinations_with_replacement(range(1,10),4) if 4 in t]
print("n",n,"|S|",len(S))
for w in words:
    os_=offs(w)
    if len(os_)<4: continue
    ins=[o for o in os_ if fm._startOf(tuple(batch._rowFor(n,w,o)),R) in S]
    if not ins or len(ins)==len(os_): continue
    for o in ins: print(w,"placements",len(os_),"o*",o,"L",o,"g",n-(o+4+int(w[-1])))
