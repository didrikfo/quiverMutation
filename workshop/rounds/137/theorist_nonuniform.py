import sys, warnings; warnings.simplefilter('ignore')
sys.path.insert(0,'.'); sys.path.insert(0,'workshop/rounds/129')
import theorist_jsilt as J
from quivermutation import mutation, procedure, arrowPaths as ap
E=[(1,2),(2,3),(3,5),(4,3)]
for name,R in (('joint abd+cd',[{(0,1,2):1,(3,2):1}]),('split',[{(0,1,2):1},{(3,2):1}])):
    A=J.build(E,R)
    print(name, 'inj',mutation.isTiltingAtVertex(A,3),'gate',mutation.mutationIsPossibleAtVertex(A,3),
          'illegal',[ap.isIllegalRelation(A.quiver,r) for r in procedure.relationsFrom(A)])
