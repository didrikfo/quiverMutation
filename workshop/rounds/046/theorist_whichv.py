"""For n=7: which (v, direction) give an LNA, and how many events are admissible per vertex."""
import collections
from quivermutation import lnaMoves as lm, mutation, nakayama, pathAlgebra
import importlib.util, sys
n = 7
spec = importlib.util.spec_from_file_location("x", "workshop/rounds/046/theorist_locality.py")
def allLNAs(n):
    out=[]
    def rec(j,cur):
        if j==n-2:
            if lm.isAdmissible(n,cur): out.append(list(cur))
            return
        for a in [0]+list(range(2,n-1-j+1)):
            cur.append(a); rels=lm.relationsOf(cur)
            if all(e[0]+e[1]<l[0]+l[1] for e,l in zip(rels,rels[1:])): rec(j+1,cur)
            cur.pop()
    rec(0,[]); return out
adm=collections.Counter(); lna=collections.Counter()
for rl in allLNAs(n):
    alg=nakayama.LinearNakayamaAlgebra(n,rl); dual=lm._quiet(pathAlgebra.dualPathAlgebra,alg)
    for v in range(1,n+1):
        for s in (v,-v):
            if not lm._quiet(mutation.mutationIsPossibleAtVertex,alg if s>0 else dual,v): continue
            adm[s]+=1
            nx_=lm._quiet(mutation.quiverMutationAtVertices,lm._copy(alg),[s])
            if nx_ is not None and lm.asRelLengths(nx_,n) is not None: lna[s]+=1
print("admissible",sorted(adm.items())); print("LNA result",sorted(lna.items()))
