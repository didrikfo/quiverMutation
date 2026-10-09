"""T5: socle reading of AI 2.32(b).  J_i = ker(x -> (x b)_b) = Hom(S_v, e_iA) = {x in e_iAe_v : x a in I for ALL arrows a}.
For two parents (E-066 step 7 at v=4; E-078 'long' square at n=5, v=d) compute, for every i != v:
  dimV  = dim e_iAe_v          (paths i~>v mod I)
  J_map = dimV - rank of p -> (p b)_b over the arrows b out of v   (tiltingPlus map)
  J_soc = dim {x : x a = 0 mod I for every arrow a of the quiver}  (socle copies of S_v, own matrix, all arrows)
  eul   = C[i,v] - sum_b C[i,t(b)]   (Cartan numbers) ; coker = J - eul ; cokerDirect = sum_b dim e_iAe_t(b) - rank
and the Cartan defect of the repo's rewrite vs R C R^T (row v off-diagonal should be -J_i, E-095).
Also checks the relations are in rad^2 (so rad P_v / rad^2 has one top per out-arrow)."""
import sys; sys.path.insert(0,'.'); sys.argv=['x']
src=open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n","\n")
exec(compile(src,'h015','exec'))
sys.path.insert(0,'workshop/rounds/011')
def rankM(rows, cols):
    import sympy
    if not rows: return 0
    M=sympy.Matrix([[r.get(c,0) for c in cols] for r in rows]); return M.rank()
def analyse(alg, v, label):
    q=alg.quiver; rels=procedure.relationsFrom(alg); verts=sorted(alg.vertices())
    outs=ap.arrowsOutOf(q,v); allarr=sorted(q.edges(keys=True))
    print("==",label,"v=",v,"out-arrows",outs,"rel lengths",[len(r) if not isinstance(r,dict) else None for r in rels][:6])
    piv=lambda i,j: ap.idealBasis(q,rels,i,j)
    C=cartan(alg); idx={x:k for k,x in enumerate(verts)}
    def cart(i,j): return C[idx[j],idx[i]] if False else C[idx[i],idx[j]]
    R=rplus(alg,v,verts)
    ch=reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg,v)); Cc=cartan(ch)
    ok=True
    for i in verts:
        if i==v: continue
        P=ap.allPathsBetween(q,i,v)
        if not P: continue
        dimV=len(P)-len(piv(i,v))
        def rowsFor(arrs):
            rows=[]
            for p in P:
                row={}
                for b in arrs:
                    for kk,val in ap.reduceAgainstPivots(ap.combination([p+(b,)]),piv(i,b[1])).items(): row[(b,kk)]=val
                rows.append(row)
            return rows
        cols=lambda rows: sorted({c for r in rows for c in r},key=str)
        r1=rowsFor(outs); r2=rowsFor([a for a in allarr if a[0]==v])  # arrows with other tails give x a = 0 identically
        Jmap=dimV-rank(r1); Jsoc=dimV-rank(r2)
        tgt=sum(len(ap.allPathsBetween(q,i,b[1]))-len(piv(i,b[1])) for b in outs)
        coker=tgt-rank(r1)
        print(" i=%s dimV=%d J_map=%d J_soc=%d dimTarget=%d coker=%d  eulerCheck %d == %d"%(i,dimV,Jmap,Jsoc,tgt,coker,dimV-tgt,Jmap-coker))
        assert Jmap==Jsoc and dimV-tgt==Jmap-coker
    R2=R.dot(C).dot(R.T)
    verts_=verts
    diff=[(verts_[a],verts_[b],int(Cc[a][b]-R2[a][b])) for a in range(len(verts_)) for b in range(len(verts_)) if Cc[a][b]!=R2[a][b]]
    print(" Cartan(child) - R C R^T nonzero entries:",diff)
# E-066 step 7 parent
lna=nk.LinearNakayamaAlgebra.fromClassName("03033030"); alg=lna.relationDual()
for w in [4,6,4,6,9,4]: alg=reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg,w))
print("gate at 4:",mutation.mutationIsPossibleAtVertex(alg,4),"tiltingPlus:",tiltingPlus(alg.quiver,procedure.relationsFrom(alg),4))
analyse(alg,4,"E-066 step-7 parent")
# E-078 n=5 long square
import importlib.util
A=pathAlgebra.PathAlgebra(); a,b,c,d,e=1,2,3,4,5
A.add_vertices_from([a,b,c,d,e])
for x,y in [(a,b),(a,c),(b,d),(c,d),(d,e)]: A.add_arrow(x,y)
A.add_rel([[a,b,d,e],[a,c,d,e]])
analyse(A,d,"E-078 abde=acde")
