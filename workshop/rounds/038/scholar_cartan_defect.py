"""Round 038 (scholar): for gate-admitted J != 0 hand algebras, C(mutated) - r C r^T (Ladkani Prop 3.6) vs sum dim J_i.
Usage (repo root): .venv/bin/python workshop/rounds/038/scholar_cartan_defect.py"""
import sys
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
_a = ARGV; exec(compile(src, 'bd', 'exec')); _a = ARGV
a3 = [(1,2),(1,3),(1,4),(2,5),(3,5),(4,5),(5,6)]
F = {'T1 (d,J)=(3,1)': (a3, [[[1,2,5,6],[1,3,5,6]]], 6, 5),
     'T2 (3,2)': (a3, [[[1,2,5,6],[1,3,5,6]],[[1,3,5,6],[1,4,5,6]]], 6, 5),
     'E-080 square a-b,a-c,b-d,c-d,d-e abde=acde (2,1) at d': ([(1,2),(1,3),(2,4),(3,4),(4,5)], [[[1,2,4,5],[1,3,4,5]]], 5, 4)}
for name, (arr, rl, n, v) in F.items():
    A = build(arr, rl, list(range(1, n+1)))
    ok = mutation.mutationIsPossibleAtVertex(A, v)
    B = mutation.quiverMutationAtVertex(A, v)
    CA = np.array(invariants.cartanMatrix(A), dtype=int); CB = np.array(invariants.cartanMatrix(B), dtype=int)
    verts = sorted(A.vertices()); k = verts.index(v)
    r = np.eye(n, dtype=int)
    for j, w in enumerate(verts):
        r[k, j] = -(j == k) + sum(1 for (a, b) in arr if a == v and b == w)
    print(name, 'gate', ok, 'dJ', {i: x for i, x in dims(A, v).items()} if 'dims' in globals() else '')
    for rr in (r, r.T):
        D = CB - rr @ CA @ rr.T
        print(' defect (CB - r CA r^T) nonzero entries:', {(i, j): int(D[i, j]) for i in range(n) for j in range(n) if D[i, j]})
    print(' CA\n', CA, '\n CB\n', CB)
