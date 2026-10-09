"""Round 041 (theorist): the n = 4 algebra A4: arrows 1->2,1->3,2->3,2->4,3->4, relation 1-2-3-4 = 1-3-4, v = 3.
Prints C_A, C' = r C_A r^T, J, C_B (actual child), characteristic polynomials of Phi = -C^{-1}C^T and the child's quiver/relations.
Usage: theorist_example.py (repo root)"""
import sys
import numpy as np, sympy as sp
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
A = build([(1,2),(1,3),(2,3),(2,4),(3,4)], [[[1,2,3,4],[1,3,4]]], [1,2,3,4]); v = 3
print('gate', mutation.mutationIsPossibleAtVertex(A, v), 'J', perI(A, v))
raw = mutation.quiverMutationAtVertex(A, v); ch = reduction.reducePathAlgebra(raw)
print('child arrows', sorted(ch.quiver.edges(keys=True)), 'rels', procedure.relationsFrom(ch))
CA = np.array(invariants.cartanMatrix(A), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
print('C_A\n', CA, '\nC_B\n', CB)
x = sp.symbols('x')
for nm, C in (('A', CA), ('B', CB)):
    M = sp.Matrix(C.tolist()); print(nm, 'det C', M.det(), 'chi =', sp.factor(sp.expand((x*M + M.T).det()/M.det())))
print('keys equal', search._coxeterKeyOrNone(A) == search._coxeterKeyOrNone(ch), search._coxeterKeyOrNone(A))
