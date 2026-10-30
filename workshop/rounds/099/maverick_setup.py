"""Round 099 maverick: structure of the n = 12 140-pair: Phi, S, commutant lattice, Jordan form at -1, roots of the degree-8 factor.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/099/maverick_setup.py"""
import sys; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
N = 12
a = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); b = (2, 2, 3, 0, 3, 0, 2, 2, 3, 0)
Ca, Cb = cartan(N, a), cartan(N, b)
Pa, Pb = phi(Ca), phi(Cb)
Sa, Sb = Ca + Ca.T, Cb + Cb.T
print("sig Sa", sig(Sa), "Sb", sig(Sb), "det", Sa.det(), Sb.det())
for P in (Pa, Pb):
    print("rank(Phi+1)", (P + sp.eye(N)).rank(), "rank(Phi^2+1)", (P**2 + sp.eye(N)).rank(), "rank (Phi+1)^2", ((P + sp.eye(N))**2).rank())
x = sp.symbols('x')
g = x**8 - x**7 - 2*x**6 - x**5 + x**4 - x**3 - 2*x**2 - x + 1
print("roots of g:", [complex(r) for r in sp.Poly(g, x).nroots()])
print("Phi order-ish: Phi^k=I?", [k for k in (1,2,3,4,5,6,8,10,12,20,24,30) if (Pa**k).is_Identity or (Pa**k) == sp.eye(N)])
