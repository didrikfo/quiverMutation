"""Test the local identity R: (a,b,b,d) -> (a-1,b,d+1) (four consecutive relation starts, d >= b, 2 <= a <= b, b >= 3)
as a one-step double-mutation neighbour at every interior placement, n = 14; also the 3-run form (a,b,b) -> (a-1,b) (nothing after).
  timeout 10m .venv/bin/python workshop/rounds/017/theorist_rrule.py [N]
"""
import sys
sys.path.insert(0,'.')
from quivermutation import freeMoves as fm, doubleMutation, lnaMoves
n=int(sys.argv[1]) if len(sys.argv)>1 else 14
def mk(word,o):
    if o+len(word)>n-2: return None
    r=[0]*(n-2)
    for i,x in enumerate(word): r[o+i]=x
    return tuple(r) if lnaMoves.isAdmissible(n,r) else None
def nbrs(r): return {fm.stripLengthTwo(tuple(v)) for v,_ in doubleMutation.rewritesOf(n,r)}
tot=ok=0; bad=[]
for b in range(3,8):
  for a in range(2,b+1):
    for d in range(b,10):
      for o in range(0,n-2):
        r=mk((a,b,b,d),o)
        if r is None: continue
        tgt=mk((a-1,b,d+1),o)  # first relation keeps its start; shifted by (d+1) at o+2? placed below
        # predicted row: starts o (a-1), o+1 (b), o+2 -> start moves to o+1? the d-relation now starts at o+1+... use intervals instead
        iv=sorted([(o+1,o+1+a),(o+2,o+2+b),(o+3,o+3+b),(o+4,o+4+d)])  # 1-based start vertices
        new=[(o+1,o+a),(o+2,o+2+b),(o+3,o+4+d)]
        pred=doubleMutation.relLengthsOf(n,new)
        tot+=1
        if pred is not None and fm.stripLengthTwo(pred) in nbrs(r): ok+=1
        else: bad.append((a,b,d,o,pred is None))
print("n",n,"cases",tot,"predicted neighbour present",ok,"missing",len(bad)); from collections import Counter
print("missing by (a,b,d), pred-is-None flag:",sorted(Counter((x[:3],x[4]) for x in bad).items()))
# three-run form (a,b,b) -> (a-1,b), 3 <= a <= b
tot=ok=0; bad=[]
for b in range(3,9):
  for a in range(3,b+1):
    for o in range(0,n-2):
        r=mk((a,b,b),o)
        if r is None: continue
        pred=doubleMutation.relLengthsOf(n,[(o+1,o+a),(o+2,o+2+b)]); tot+=1
        if pred is not None and fm.stripLengthTwo(pred) in nbrs(r): ok+=1
        else: bad.append((a,b,o))
print("three-run: cases",tot,"present",ok,"missing",bad[:10])
