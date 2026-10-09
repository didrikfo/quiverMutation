"""Round 015 theorist: trace step 7 of the mutation at vertex 3 of the first n = 8 class 2 'M' parent."""
import sys
exec(open('workshop/rounds/015/theorist_cartan.py').read().split("for kind, t in lines:")[0].replace("n, c = int(sys.argv[1]), int(sys.argv[2])", "n, c = 8, 2").replace("sys.argv = ['x']", "sys.argv=['x']"))
from quivermutation import procedure as P
kind, t = lines[1]
d, rels, v, path = t
alg = classes[base][path[0]]
for w in path[1:]:
    alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
rl = P.relationsFrom(alg)
print('relations (combos):'); [print('  ', r) for r in rl]
q, nr = P.mutateAtVertex(alg.quiver, rl, v)
print('raw mutated relations:'); [print('  ', r) for r in nr]
outRel = [r for r in rl if ap.source(r) == v]
print('out relations of v:', outRel)
# does the combination 3*>2>7 + 3*>6>7 belong to the ideal of the true End algebra? compute kernel directly
# --- trace step 7
outArrows = ap.arrowsOutOf(alg.quiver, v)
inArrows = ap.arrowsInto(alg.quiver, v)
print('out arrows', outArrows, 'in arrows', inArrows)
orig = P._kernelOverIdeal
def traced(quiver, relations, candidates, outArrows, shadow):
    res = orig(quiver, relations, candidates, outArrows, shadow)
    print('TARGET candidates', candidates); print('  shadow', shadow); print('  kernel', res)
    return res
P._kernelOverIdeal = traced
q, nr = P.mutateAtVertex(alg.quiver, rl, v)
print('--- residues at target 7')
for el in ({((1, 2, 0), (2, 7, 0)): 1}, {((1, 5, 0), (5, 6, 0), (6, 7, 0)): -1}):
    print(el, '->', P._reduceAgainstIdeal(alg.quiver, rl, el))
print('ideal basis 1->7', ap.idealBasis(alg.quiver, rl, 1, 7))
print('all paths 1->7', ap.allPathsBetween(alg.quiver, 1, 7))
