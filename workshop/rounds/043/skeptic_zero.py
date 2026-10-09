"""Round 043 (skeptic): re-build, from the pickle written by skeptic_c2.py (the first D == 0 steps), the parent and child by hand
(arrows + relations only) and run: gate, tiltingPlus (Ladkani 2.3(c)), Cartan congruence C_B = R C_A R^T (Ladkani 3.6), Coxeter keys,
and J_i.   Usage: skeptic_zero.py /path/zero_n_c.pkl"""
import sys, pickle, ast
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
def mk(edges, reprr, nodes):
    arrows = [(a, b) for a, b, k in edges]; rels = []
    for d in ast.literal_eval(reprr):
        ps = [[p[0][0]] + [a[1] for a in p] for p in d]
        rels.append(ps)
    return build(arrows, rels, nodes)
for z in pickle.load(open(ARGV[1], 'rb')):
    N = sorted({x for e in z['pe'] for x in e[:2]}); v = z['v']
    A = mk(z['pe'], z['pr'], N); B = mk(z['ce'], z['cr'], N)
    rels = procedure.relationsFrom(A)
    print('parent arrows', [e[:2] for e in z['pe']], 'rels', z['pr'][:200])
    print(' v', v, 'J', z['J'], 'gate', mutation.mutationIsPossibleAtVertex(A, v), 'tiltingPlus', bool(tiltingPlus(A.quiver, rels, v)))
    raw = mutation.quiverMutationAtVertex(A, v); ch = reduction.reducePathAlgebra(raw)
    print(' key(A)', search._coxeterKeyOrNone(A), 'key(child)', search._coxeterKeyOrNone(ch), 'hand child key', search._coxeterKeyOrNone(B),
          'canonical same', fingerprint.canonicalKey(ch) == fingerprint.canonicalKey(B))
    R = rplus(A, v, N); CA = cartan(A); CB = cartan(ch)
    print(' Cartan congruent (R C_A R^T == C_B):', bool((R.dot(CA).dot(R.T) == CB).all()), ' dim A', int(CA.sum()), 'dim B', int(CB.sum()))
