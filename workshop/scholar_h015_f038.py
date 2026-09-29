import sys; sys.path.insert(0,'.'); sys.argv=['x']
import importlib.util
src=open('workshop/scholar_h015.py').read().replace("\nmain()\n","\n")
exec(compile(src,'h015','exec'))
lna=nk.LinearNakayamaAlgebra.fromClassName("03033030")
import os; alg=lna.relationDual() if os.environ.get("RD") else pathAlgebra.dualPathAlgebra(lna)
print(len(alg.vertices()))
for step,v in enumerate([4,6,4,6,9,4,4,6],1):
    ok=mutation.mutationIsPossibleAtVertex(alg,v)
    rels=procedure.relationsFrom(alg)
    tilt=tiltingPlus(alg.quiver,rels,v) if ok else None
    child=reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg,v)) if ok else None
    verts=sorted(alg.vertices()); R=rplus(alg,v,verts)
    cong=bool((R.dot(cartan(alg)).dot(R.T)==cartan(child)).all()) if ok else None
    print(step,v,'gate',ok,'tilt',tilt,'cong',cong,'key same',search._coxeterKeyOrNone(child)==search._coxeterKeyOrNone(alg) if ok else None)
    alg=child
