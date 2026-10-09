"""Round 051: replay the five F-037 depth-7 paths 34504030 -> 50505000 (n = 10) edge by edge:
gate, tiltingPlus (J = 0), Coxeter key kept.  Reuses the replay of toolsmith_witness.py (047).
Usage: .venv/bin/python workshop/rounds/051/experimentalist_f037replay.py"""
import sys, copy
sys.path.insert(0, '.')
from quivermutation import lnaMoves as lm, nakayama as nk, mutation, procedure, reduction, search, pathAlgebra
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'))
paths = [[4,2,1,1,2,2,4],[4,1,2,2,2,4,1],[4,1,2,2,2,1,4],[4,1,2,2,1,2,4],[4,1,1,2,2,2,4]]
def replay(alg0, path):
    alg = copy.deepcopy(alg0); base = search._coxeterKeyOrNone(alg); edges = []
    for v in path:
        gate = mutation.mutationIsPossibleAtVertex(alg, v)
        J0 = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v)) if gate else None
        if not gate: edges.append((v, False, None, None)); break
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(alg), v))
        edges.append((v, gate, J0, search._coxeterKeyOrNone(alg) in (None, base)))
    return alg, edges
A = nk.LinearNakayamaAlgebra(10, [3,4,5,0,4,0,3,0])
B = nk.LinearNakayamaAlgebra(10, [5,0,5,0,5,0,0,0])
tot = [0, 0, 0, 0]
for p in paths:
    alg, es = replay(A, p)
    row = lm.asRelLengths(lm._copy(alg), 10); arrive = row
    print(p, 'gate', sum(bool(e[1]) for e in es), 'J0', sum(bool(e[2]) for e in es), 'key', sum(bool(e[3]) for e in es),
          'of', len(es), 'ends on LNA row:', arrive, flush=True)
    for e in es: print('   ', e)
