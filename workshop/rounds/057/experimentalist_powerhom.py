"""Round 057: Hom(T,T[+-1]) = 0 (skeptic_tilt.stepTest, independent of tiltingPlus) on every gate+key edge of the depth-D walk tree from rows A, B (right and dual steps),
compared with tiltingPlus (J = 0).  Usage (repo root): .venv/bin/python workshop/rounds/057/experimentalist_powerhom.py N ROWA ROWB DEPTH"""
import sys, copy, collections
sys.path.insert(0, '.')
import numpy as np
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import nakayama as nk, mutation, procedure, reduction, search, pathAlgebra, arrowPaths
import importlib.util
spec = importlib.util.spec_from_file_location('am', 'workshop/rounds/054/experimentalist_amerge.py'); am = importlib.util.module_from_spec(spec); spec.loader.exec_module(am)
TP = am.loadTP()
n = int(sys.argv[1]); D = int(sys.argv[4]); tally = collections.Counter(); seen = set()
def go(a, d, base, tag):
    if d == 0 or search.hasOrientedCycle(a): return
    k = search.quiverKey(a)
    if k is not None:
        if (k, d) in seen: return
        seen.add((k, d))
    for v in sorted(a.quiver.nodes):
        if not mutation.mutationIsPossibleAtVertex(a, v): continue
        j0 = bool(TP(a.quiver, procedure.relationsFrom(a), v))
        m = mutation.quiverMutationAtVertex(copy.deepcopy(a), v)
        if any(arrowPaths.isIllegalRelation(m.quiver, r) for r in procedure.relationsFrom(m)): continue
        m = reduction.reducePathAlgebra(m)
        mk = search._coxeterKeyOrNone(m)
        if mk is not None and mk != base: continue
        Vs, res = stepTest(a, v, -1)
        tilt = not np.abs(np.array(res[-1])).sum() and not np.abs(np.array(res[1])).sum()
        tally[(tag, 'J0' if j0 else 'J!=0', 'tilt' if tilt else 'NOT tilt')] += 1
        go(m, d - 1, base, tag)
for row in sys.argv[2:4]:
    alg = nk.LinearNakayamaAlgebra(n, [int(c) for c in row]); base = search._coxeterKeyOrNone(alg)
    for dual in (False, True):
        go(pathAlgebra.dualPathAlgebra(alg) if dual else alg, D, base, row + ('d' if dual else ''))
for k, v in sorted(tally.items()): print(k, v)
print('total edges', sum(tally.values()), '| J0&tilt', sum(v for k, v in tally.items() if k[1] == 'J0' and k[2] == 'tilt'), '| J!=0', sum(v for k, v in tally.items() if k[1] == 'J!=0'))
