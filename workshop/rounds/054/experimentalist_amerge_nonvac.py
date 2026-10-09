"""Round 054: non-vacuity of the J test. From the n=10 start and its depth-1,2 forward neighbours, count
vertices where the gate admits a mutation and tiltingPlus gives J != 0 (so the J test can fail on this path).
Usage (repo root): timeout 10m .venv/bin/python workshop/rounds/054/experimentalist_amerge_nonvac.py"""
import sys, copy, importlib.util
sys.path.insert(0, '.')
spec = importlib.util.spec_from_file_location('am', 'workshop/rounds/054/experimentalist_amerge.py'); am = importlib.util.module_from_spec(spec); spec.loader.exec_module(am)
from quivermutation import nakayama as nk, mutation, procedure, pathAlgebra, reduction
tp = am.loadTP()
def scan(alg, tag):
    g = [v for v in range(1, 11) if mutation.mutationIsPossibleAtVertex(alg, v)]
    j = [v for v in g if not tp(alg.quiver, procedure.relationsFrom(alg), v)]
    print(tag, 'gate-admitted', g, 'of which J != 0:', j); return g, j
alg = nk.LinearNakayamaAlgebra(10, list(am.START))
tot = bad = 0
for tag, a in (('start', alg), ('start-dual', pathAlgebra.dualPathAlgebra(alg))):
    g, j = scan(a, tag); tot += len(g); bad += len(j)
    for v in g:
        b = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(a), v))
        g2, j2 = scan(b, '  after %s%d' % ('-' if tag.endswith('dual') else '', v)); tot += len(g2); bad += len(j2)
print('gate-admitted edges examined', tot, 'with J != 0:', bad)
# wider: start vertices of all distinct checkpoint members (n=10 LNAs), both orientations
import json
mem = sorted({tuple(json.loads(l)['member']) for l in open(am.CK)})
n = nb = 0; ex = None
for m in mem:
    a0 = nk.LinearNakayamaAlgebra(10, list(m))
    for a in (a0, pathAlgebra.dualPathAlgebra(a0)):
        for v in range(1, 11):
            if mutation.mutationIsPossibleAtVertex(a, v):
                n += 1
                if not tp(a.quiver, procedure.relationsFrom(a), v):
                    nb += 1; ex = ex or (m, a is not a0, v)
print('checkpoint members', len(mem), 'gate-admitted (member, orientation, vertex) triples', n, 'with J != 0:', nb, 'first example', ex)
