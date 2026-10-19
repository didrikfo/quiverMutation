"""Round 074 (skeptic): hand-check of ONE child of the 70 (a 'new' iso class, class 2, v = 6).  Prints path bases, Cartan matrices, the independent Hom(T,T[m]) table
(rounds/050/skeptic_tilt.py), my own Coxeter polynomials and the structure of the child.  Run from repo root: needs /tmp/claude-0/sk74/c2_s*.pkl (skeptic_failkids.py)."""
import sys, pickle
sys.path.insert(0, '.')
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from sympy import Matrix, symbols, factor
import numpy as np
D = sys.argv[1] if len(sys.argv) > 1 else '/tmp/claude-0/sk74'
want = sorted([(1, 7, 0), (3, 1, 0), (3, 4, 0), (4, 5, 0), (5, 6, 0), (6, 2, 0), (7, 6, 0)])
rec = None
for s in range(4):
    for (d, pk, v, pa, ch) in pickle.load(open('%s/c2_s%d.pkl' % (D, s), 'rb'))[1]:
        if v == 6 and sorted(pa.quiver.edges(keys=True)) == want: rec = (d, v, pa, ch)
d, v, pa, ch = rec
print('parent depth', d, 'v', v)
for name, a in (('parent', pa), ('child', ch)):
    A = Alg(a); V = sorted(a.quiver.nodes)
    print(name, 'arrows', sorted(a.quiver.edges(keys=True)))
    C = np.array([[len(A.basis(i, j)) for j in V] for i in V])      # C[i][j] = dim e_i A e_j (paths i~>j mod I; i = j gives the trivial path)
    print(name, 'dim', C.sum(), 'Cartan (row i, col j = #paths i~>j mod I):'); print(C)
    x = symbols('x'); M = Matrix(C.tolist())
    # Coxeter matrix Phi = -C^{-T} C ; poly det(xI - Phi)
    Phi = -(M.T.inv()) * M
    print(name, 'Coxeter polynomial', factor((x * Matrix.eye(len(V)) - Phi).det()))
V, res = stepTest(pa, v, -1)
Cc = np.array([[len(Alg(ch).basis(i, j)) for j in V] for i in V])
print('V', V)
for m in (-1, 0, 1):
    print('dim Hom(T_i,T_j[%d]) (row i, col j):' % m); print(np.array(res[m]))
print('End(T) Cartan (m = 0 table) equals the child Cartan:', bool((np.array(res[0]) == Cc).all()))
print('Hom(T,T[-1]) nonzero entries (i, j, dim):', [(V[i], V[j], res[-1][i][j]) for i in range(len(V)) for j in range(len(V)) if res[-1][i][j]])
print('Hom(T,T[+1]) nonzero entries:', [(V[i], V[j], res[1][i][j]) for i in range(len(V)) for j in range(len(V)) if res[1][i][j]])
