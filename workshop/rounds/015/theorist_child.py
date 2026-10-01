"""Round 015 theorist: print parent and child (arrows, relations, Cartan columns) for the i-th replayed line of n = 8 class 2.
Usage: theorist_child.py i"""
import sys
i = int(sys.argv[1])
exec(open('workshop/rounds/015/theorist_cartan.py').read().split("for kind, t in lines:")[0].replace("n, c = int(sys.argv[1]), int(sys.argv[2])", "n, c = 8, 2").replace("sys.argv = ['x']", "sys.argv=['x']"))
kind, t = lines[i]
if kind == 'R': d, rels, v, g, path = t
else: d, rels, v, path = t
alg = classes[base][path[0]]
for w in path[1:]:
    alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
raw = mutation.quiverMutationAtVertex(alg, v); ch = reduction.reducePathAlgebra(raw)
for name, a in (('parent', alg), ('child raw', raw), ('child reduced', ch)):
    print(name, 'vertices', sorted(a.vertices()), 'arrows', sorted((x, y) for x, y, *_ in a.quiver.edges(keys=True)))
    print('  rels', a.rels)
verts = sorted(alg.vertices())
print('Cartan parent'); print(cartan(alg)); print('predicted'); print(rplus(alg, v, verts).dot(cartan(alg)).dot(rplus(alg, v, verts).T)); print('child'); print(cartan(ch))
