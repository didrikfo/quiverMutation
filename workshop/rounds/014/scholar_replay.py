"""Replay the first recorded rejection of a scholar_walk.py class from its start, step by step.
Usage: scholar_replay.py n class   (reads workshop/rounds/014/scholar_walk_n{n}_c{class}.txt)
Checks independently of the BFS: every prefix step gate-admitted with child key = base key (guarded walk),
then at the last vertex gate True, tiltingPlus False, child key != base key."""
import sys, re, ast
n, c = int(sys.argv[1]), int(sys.argv[2])
line = [l for l in open('workshop/rounds/014/scholar_walk_n%d_c%d.txt' % (n, c)) if l.startswith('(') and 'False, (' in l][0]
d, rels, v, guard, path = ast.literal_eval(line.strip())
sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[c]; alg = classes[base][path[0]]
print('start', alg.rels, 'path', path[1:], 'reject vertex', v)
for step, w in enumerate(path[1:], 1):
    assert mutation.mutationIsPossibleAtVertex(alg, w)
    child = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
    assert search._coxeterKeyOrNone(child) == base, ('guard', step)
    alg = child
rl = procedure.relationsFrom(alg)
g = mutation.mutationIsPossibleAtVertex(alg, v); t = tiltingPlus(alg.quiver, rl, v)
ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, v))
print('parent: vertices', len(alg.vertices()), 'arrows', sorted((a, b) for a, b, *_ in alg.quiver.edges(keys=True)) if hasattr(alg.quiver, 'edges') else '')
print('relations', alg.rels)
print('depth', len(path) - 1, 'gate', g, 'tiltingPlus', t, 'child key == base', search._coxeterKeyOrNone(ch) == base)
print('child key', search._coxeterKeyOrNone(ch), 'base', base)
