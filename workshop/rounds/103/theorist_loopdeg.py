"""Round 103 (theorist), T5 / H-023 / E-256 a: can the library be handed a step whose dg model has a loop of degree >= 1 at i?
Parent of E-211 (a:1->2, s:2->3, b:3->4, c:4->3; monomial relations a s, s b c, c b), i = 1 (not on the cycle) and i = 3, 4 (on the
cycle 3 <-> 4; the relation c b (c then b) is a minimal relation 4 -> 4, i.e. a degree-1 loop at 4 in Oppermann model).
Prints: gate, default call, allowCyclic=True call, for each i.  usage (repository root):
  .venv/bin/python -W ignore workshop/rounds/103/theorist_loopdeg.py"""
import sys, traceback
sys.path.insert(0, '.'); sys.path.insert(0, 'tests')
import networkx as nx
from quivermutation import procedure as pr, mutation, pathAlgebra, arrowPaths as ap
from helpers import path_algebra, arrow_set, rel_set, quiet

A = [(1, 2), (2, 3), (3, 4), (4, 3)]
# relations as vertex paths: a s = 1 2 3 ; s b c = 2 3 4 3 ; c b = 4 3 4
R = [[(1, 2, 3)], [(2, 3, 4, 3)], [(4, 3, 4)]]
alg = path_algebra(A, R)
print('parent', sorted(arrow_set(alg)), sorted(rel_set(alg)))
Q = alg.quiver; rels = pr.relationsFrom(alg)
for i in (1, 3, 4):
    below = nx.descendants(Q, i) | {i}
    print('== i =', i, '| D(i) =', sorted(below), '| D(i) cyclic:', ap.hasOrientedCycle(Q.subgraph(below)),
          '| i on a cycle:', any(i in c for c in nx.simple_cycles(nx.DiGraph(Q))))
    try: print(' gate (isMutable):', pr.isMutable(Q, rels, i))
    except Exception as e: print(' gate RAISED', type(e).__name__, str(e)[:80])
    for allow in (False, True):
        try:
            Q2, R2 = pr.mutateAtVertex(Q, rels, i, allowCyclic=allow)
            print(' allowCyclic=%s -> returns, %d arrows, %d relations' % (allow, Q2.number_of_edges(), len(R2)))
        except Exception as e:
            print(' allowCyclic=%s -> %s: %s' % (allow, type(e).__name__, str(e)[:90]))
