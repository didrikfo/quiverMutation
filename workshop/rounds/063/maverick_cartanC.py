"""Round 063: Cartan-level Candidate C on the n = 10 Phi^18 group (E-172): for one representative per class,
order of Phi, [S] = -Phi, smallest (a,b) with [S]^a = (-1)^b I, signature/det of C+C^T, gldim-free numerics.
usage (repo root): .venv/bin/python workshop/rounds/063/maverick_cartanC.py"""
import sympy
from quivermutation import coxeterTables as ct
n = 10
reps = {"class0 00000030": (0,0,0,0,0,0,3,0,0,0)[:9], "34504030": (3,4,5,0,4,0,3,0,0,0)[:9] , "50505000": (5,0,5,0,5,0,0,0,0,0)[:9], "90000000": (9,0,0,0,0,0,0,0,0,0)[:9]}
# relation-vector strings in the repo have 8 digits (n-2); take them as written
strs = {"class0 00000030": "00000030", "34504030": "34504030", "50505000": "50505000", "90000000": "90000000"}
for name, s in strs.items():
    rl = tuple(int(c) for c in s)
    C = sympy.Matrix(ct.lnaCartanMatrix(n, rl)); I = sympy.eye(n)
    Phi = -(C.inv().T) * C; S = -Phi
    Ps = Phi * 1
    ordr = next(e for e in range(1, 61) if Phi**e == I)
    pm = [a for a in range(1, 19) if Phi**a == -I]
    sym = C + C.T
    ev = sorted(float(x) for x in sym.evalf().eigenvals(multiple=True)) if False else None
    import numpy as np
    w = np.linalg.eigvalsh(np.array(sym.tolist(), dtype=float))
    print(name, "ord(Phi)=", ordr, "Phi^a=-I for a in", pm, "[S]^a=I for a in", [a for a in range(1,19) if S**a == I],
          "| C+C^T eig: neg", int((w < -1e-9).sum()), "zero", int((abs(w) < 1e-9).sum()), "pos", int((w > 1e-9).sum()), "min %.3f" % w.min(),
          "det C+C^T", sym.det())
# D_10 reference: Cartan of kD_10 (any orientation)
edges = [(i, i+1) for i in range(8)] + [(7, 9)]
C = sympy.eye(n)
for i, j in edges: C[i, j] = 1
Phi = -(C.inv().T) * C
print("kD_10", "Phi^9=-I:", Phi**9 == -sympy.eye(n), "Phi^18=I:", Phi**18 == sympy.eye(n), "char poly", sympy.factor(Phi.charpoly().as_expr()))
# F-045 check: G = C^-1 + C^-T positive definite, det, SNF (D_r has det 4); and kA_10/rad^9 = digits 90000000
from sympy.matrices.normalforms import smith_normal_form
for s in strs.values():
    C = sympy.Matrix(ct.lnaCartanMatrix(n, tuple(int(c) for c in s))); G = C.inv() + C.inv().T
    w = np.linalg.eigvalsh(np.array(G.tolist(), dtype=float))
    snf = smith_normal_form(G, domain=sympy.ZZ)
    print(s, "G neg", int((w < -1e-9).sum()), "zero", int((abs(w) < 1e-9).sum()), "det", G.det(), "SNF diag", [abs(int(snf[i, i])) for i in range(n)])
