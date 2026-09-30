"""E-032 ALARM step 7: parent data and the Aihara-Iyama 2.32(b)-type witness.

Replays [4,6,4,6,9,4,4,6] from the relation dual of 03033030 (RD=1 as in
scholar_h015_f038.py), stops at step 7 and prints the parent quiver, relations,
and every (i, path) in e_i A e_4 (paths i ~> 4 nonzero mod I) whose composite with
EVERY arrow out of 4 is zero mod I: a kernel element of the map
  Hom(P_i-side, P_4) -> prod_beta Hom(., P_head(beta)),  p |-> (p beta)_beta,
which is Hom(g, D) of AI Thm 2.32(b) with g: D_4 -> P_4 built from the arrows at 4.
"""
import sys; sys.path.insert(0,'.'); sys.argv=['x']
src=open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n","\n")
exec(compile(src,'h015','exec'))
lna=nk.LinearNakayamaAlgebra.fromClassName("03033030")
alg=lna.relationDual()
for v in [4,6,4,6,9,4]:
    alg=reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg,v))
k=4
q=alg.quiver; rels=procedure.relationsFrom(alg)
print("arrows:",sorted(q.edges(keys=True)) if hasattr(q,'edges') else q)
print("relations:",rels)
outs=ap.arrowsOutOf(q,k); print("arrows out of 4:",outs)
print("ins of 4:",[e for e in q.edges(keys=True) if e[1]==k] if hasattr(q,'edges') else None)
def piv(i,j): return ap.idealBasis(q,rels,i,j)
for i in q.nodes:
    if i==k: continue
    P=ap.allPathsBetween(q,i,k)
    if not P: continue
    nz=[p for p in P if ap.reduceAgainstPivots(ap.combination([p]),piv(i,k))]
    dimV=len(P)-len(piv(i,k))
    rows=[]
    for p in P:
        row={}
        for b in outs:
            res=ap.reduceAgainstPivots(ap.combination([p+(b,)]),piv(i,b[1]))
            for kk,v in res.items(): row[(b,kk)]=v
        rows.append(row)
    print("i",i,"paths",len(P),"dim e_iAe_4",dimV,"rank of map",rank(rows),"-> kernel" if rank(rows)<dimV else "ok")
    if rank(rows)<dimV:
        for p,r in zip(P,rows):
            if not r: print("   path killed by all arrows out of 4:",p)
