"""Round 139 (scholar, response to referee): check the explicit witness of claim (3)
directly, no search: V = [-(e_1+..+e_{n-1}), e_n, e_1-e_n, e_2, .., e_{n-2}] satisfies
V^T E_L V = E_D exactly (E = C^{-1}; L = N_n(n-1), D = kD_n as in scholar_isometry.py),
det V = -1. Reproduces the referee's n = 4..20 check (17 of 17).
usage (from repo root): .venv/bin/python -I workshop/rounds/139/scholar_shape.py [NMAX]"""
import sys; sys.path.insert(0,'.')
import sympy
from quivermutation.coxeterTables import lnaCartanMatrix
def hered(n):
    arrows=[(0,2),(1,2)]+[(i,i+1) for i in range(2,n-1)]
    M=sympy.eye(n); reach={i:{i} for i in range(n)}
    for _ in range(n):
        for s,t in arrows:
            for i in range(n):
                if s in reach[i]: reach[i].add(t)
    for s in range(n):
        for t in reach[s]: M[t,s]=1
    return M
ok=0
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
for n in range(4, nmax + 1):
    e=lambda i: sympy.Matrix([1 if j==i-1 else 0 for j in range(n)])
    cols=[-sum((e(i) for i in range(1,n)),sympy.zeros(n,1)), e(n), e(1)-e(n)]+[e(i) for i in range(2,n-1)]
    V=sympy.Matrix.hstack(*cols)
    EL=sympy.Matrix(lnaCartanMatrix(n,[n-1]+[0]*(n-3))).inv(); ED=hered(n).inv()
    good = V.T*EL*V==ED; ok+=good
    print(n, good, V.det())
print(f'witness holds at {ok} of {nmax - 3} values n = 4..{nmax}')
