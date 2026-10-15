"""Round 067 (toolsmith), T10/T5 item 1: oriented cycle / loop scan of every algebra along the 5 printed E-169 n = 10 paths
(the parent of each edge = the algebra mutated, dual for R moves; plus the end).  The n = 7 paths of E-157/E-160/E-163 were scanned in
round 054 (E-168, workshop/rounds/054/theorist_loops_out.txt): 370 edges, 0 loops, 0 oriented cycles; not rerun here.
Usage (repo root): .venv/bin/python workshop/rounds/067/toolsmith_cycscan.py"""
import sys, copy
sys.path.insert(0, '.')
from quivermutation import nakayama as nk, mutation, reduction, pathAlgebra, arrowPaths as ap
START = [0, 5, 0, 4, 0, 3, 3, 0]
PATHS = [[2, 1, 6, 5, -4, -10, -8], [-10, -8, -4, -5, 2, 1, 6], [4, 3, 2, 1, 7, 6, 9], [4, 3, 7, 2, 1, 6, 9], [4, 3, 7, 9, 2, 1, 6]]
tot = dict(paths=0, parents=0, cyclic=0, ends=0, endcyclic=0, parallel=0)
for p in PATHS:
    alg = nk.LinearNakayamaAlgebra(10, list(START)); tot['paths'] += 1
    for s in p:
        v = abs(s); a = pathAlgebra.dualPathAlgebra(alg) if s < 0 else alg
        tot['parents'] += 1; tot['cyclic'] += ap.hasOrientedCycle(a.quiver)
        tot['parallel'] += any(a.quiver.number_of_edges(x, y) > 1 for x, y in a.quiver.edges())
        a2 = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(a), v))
        alg = pathAlgebra.dualPathAlgebra(a2) if s < 0 else a2
    tot['ends'] += 1; tot['endcyclic'] += ap.hasOrientedCycle(alg.quiver)
print(tot)
