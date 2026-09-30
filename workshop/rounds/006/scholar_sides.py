"""Both sides of the AI 2.32 / Ladkani test along the E-032 ALARM path (RD=1 start).
right: p: i~>k composed with each arrow out of k (tiltingPlus, the repo's test);
left : p: k~>j preceded by each arrow into k (mirror test; AI 2.32(a)-type).
Prints per step: gate, right, left, key same. Seconds."""
import sys; sys.path.insert(0,'.'); sys.argv=['x']
src=open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n","\n")
exec(compile(src,'h015','exec'))
def tiltingLeft(q,rels,k):
    ins=[e for e in q.edges(keys=True) if e[1]==k]
    piv=lambda i,j: ap.idealBasis(q,rels,i,j)
    for j in q.nodes:
        if j==k: continue
        P=ap.allPathsBetween(q,k,j)
        if not P: continue
        dimV=len(P)-len(piv(k,j)); rows=[]
        for p in P:
            row={}
            for b in ins:
                res=ap.reduceAgainstPivots(ap.combination([(b,)+p]),piv(b[0],j))
                for kk,v in res.items(): row[(b,kk)]=v
            rows.append(row)
        if rank(rows)<dimV: return False
    return True
alg=nk.LinearNakayamaAlgebra.fromClassName("03033030").relationDual()
for step,v in enumerate([4,6,4,6,9,4,4,6],1):
    ok=mutation.mutationIsPossibleAtVertex(alg,v); rels=procedure.relationsFrom(alg)
    child=reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg,v))
    print(step,v,'gate',ok,'right',tiltingPlus(alg.quiver,rels,v),'left',tiltingLeft(alg.quiver,rels,v),
          'key same',search._coxeterKeyOrNone(child)==search._coxeterKeyOrNone(alg))
    alg=child
