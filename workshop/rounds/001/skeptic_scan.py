import sys,itertools
sys.path.insert(0,'.')
import batch
from quivermutation import coxeterTables as ct
def sums(n,core):
    keys={}
    for o in range(n):
        r=batch._rowFor(n,core,o)
        if r: keys.setdefault(ct.lnaCoxeterKey(n,tuple(r)),[]).append(o)
    offs=sum(1 for o in range(n) if batch._rowFor(n,core,o))
    gs=[g for g in keys.values() if len(g)>1]
    bad=any(len(g)>2 for g in gs)
    return sorted({sum(g) for g in gs}),bad,offs
n=int(sys.argv[1])
out=[]
for L in range(1,5):
  for w in itertools.product('03456789',repeat=L):
    c=''.join(w)
    if c[0]=='0' or c[-1]=='0' or '00' in c: continue
    if not batch._rowFor(n,c,0): continue
    s,bad,offs=sums(n,c)
    out.append((c,[n-x for x in s],bad,offs))
for c,k,b,o in out:
    if len(k)!=1 or b: print(c,'n-sums',k,'triple' if b else '',o)
print(len(out),'cores;',sum(1 for c,k,b,o in out if len(k)==1 and not b),'with exactly one centre;',[c for c,k,b,o in out if not k])
