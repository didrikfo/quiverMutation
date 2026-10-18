"""Round 073 (toolsmith): do Cartan-congruence invariants beyond the Coxeter polynomial separate the certified-inequivalent F-010 pair (n = 9)?
Pair 3060000 vs 3304000 (LNA class names, relation lengths); invariants of C (Cartan, exact): SNF(C), SNF(C+C^T), SNF(C-C^T), det, signature of C+C^T,
charpoly of Coxeter Phi = -C^-T C. All are invariants of Z-congruence C -> P^T C P (P in GL_n(Z)), hence of derived equivalence.
  .venv/bin/python -W ignore workshop/rounds/073/toolsmith_f010inv.py"""
import sys; sys.path.insert(0, '.')
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import nakayama as nk, invariants
def inv(name):
    a = nk.LinearNakayamaAlgebra.fromClassName(name)
    C = sp.Matrix(invariants.cartanMatrix(a, exact=True).tolist())
    snf = lambda M: tuple(sorted(abs(int(smith_normal_form(M, domain=sp.ZZ)[i, i])) for i in range(M.rows)))
    S = C + C.T
    ev = [complex(e).real for e in S.evalf().eigenvals(multiple=True)]
    sig = (sum(e > 1e-9 for e in ev), sum(e < -1e-9 for e in ev))
    phi = (-C.T.inv() * C) if False else (-(C.T).inv() * C)
    return dict(det=C.det(), snfC=snf(C), snfSym=snf(S), snfAnti=snf(C - C.T), sig=sig, cox=sp.factor(phi.charpoly().as_expr()))
for nm in sys.argv[1:] or ['3060000', '3304000']:
    try: print(nm, inv(nm))
    except Exception as e: print(nm, 'ERR', repr(e))
