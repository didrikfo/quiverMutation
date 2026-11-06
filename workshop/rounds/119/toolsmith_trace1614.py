"""Round 119 (toolsmith), H-030: trace seed 1614 (N, n = 5) in the repo's isMutable / mutateAtVertex / reduce.
usage: toolsmith_trace1614.py   (run from repo root)"""
import sys, random
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
import skeptic_conv3 as c3
from skeptic_conv3 import AlgG, repo_input, repo_children
import skeptic_convention as sc

A = c3.gen('N', 5, random.Random(1614))
print('arrows', [(a, b) for a, b, _ in A.arrows]); print('rels', A.rels, 'minimal', c3.is_minimal(A))
q, rels = repo_input(A)
print('repo rels', rels)
for i in range(5):
    invR = sc.invariants(A, i, 'R'); invL = sc.invariants(A, i, 'L')
    kids, gF, gR = repo_children(A, i)
    print('v', i, 'tilting F', invR[0], 'R', invL[0], 'gate F', gF, 'R', gR)

# ---- the gate, step by step (F, vertex idx 3 = repo vertex 4): gate False, tilting True ----
from quivermutation import procedure, arrowPaths as ap
def gate_trace(q, rels, v):
    out = ap.arrowsOutOf(q, v)
    print('gate at', v, 'out arrows', out)
    for s in q.nodes:
        if s == v: continue
        for path in ap.allPathsBetween(q, s, v):
            z = ap.isInIdeal(q, rels, ap.combination([path]))
            ext = [(a, ap.isInIdeal(q, rels, ap.combination([path + (a,)]))) for a in out]
            print('  path', path, 'zero' if z else 'nonzero', 'extensions zero?', [e[1] for e in ext])
gate_trace(q, rels, 4)
print('--- same, relation 2 deleted')
gate_trace(q, rels[:1], 4)
