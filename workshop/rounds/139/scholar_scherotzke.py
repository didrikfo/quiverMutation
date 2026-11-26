"""Round 139 (scholar): a K0-level check of the shift statement in Scherotzke,
arXiv:0906.4987v1, Theorem 5.7 (source l. 1145): for A = kA_n/(path of length n-1),
n >= 4, the AR quiver of D^b(A) is Z[D_n], [-1] acts as the identity on D_n for n
even and as the involution for n odd, and A ~ kD_n.
On K0, tau acts by the Coxeter matrix Phi and [1] by -1. In D^b(kD_n), tau^{-(n-1)}
= [1] o sigma (h = 2n-2, sigma induced by -w0: trivial iff n even). So the theorem
predicts: Phi^(2n-2) = I, and rank(Phi^(n-1) + I) = 0 for n even, 1 for n odd
(sigma swaps two classes). This is a NECESSARY test only (K0, not the category).
Controls, where the prediction must FAIL: A_n-type control N_n(2) (rad^2 = 0,
~der kA_n, h = n+1), and E_6, E_7,
E_8 as N_6(4), N_7(5), N_8(6) (family N_{r+2}(r) ~ [2,3,r-1], LMR 2112.15587 Thm 4.4)
and N_7(3) (E_7).
usage (from repo root): .venv/bin/python -I workshop/rounds/139/scholar_scherotzke.py [NMAX]"""
import sys
sys.path.insert(0, '.')
import sympy
from quivermutation.coxeterTables import lnaCartanMatrix


def rels(n, r):
    return [r] * (n - r) + [0] * (r - 2)


def predicted(n, r):
    C = sympy.Matrix(lnaCartanMatrix(n, rels(n, r)))
    Phi = -C.T * C.inv()
    I = sympy.eye(n)
    half = Phi ** (n - 1)
    ok_full = half * half == I
    rk = (half + I).rank()
    return ok_full and rk == (0 if n % 2 == 0 else 1), ok_full, rk


nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
fails = 0
for n in range(4, nmax + 1):
    ok, full, rk = predicted(n, n - 1)
    fails += not ok
    print(f'N_{n}({n-1}): Phi^(2n-2)=I {full}, rank(Phi^(n-1)+I)={rk}  {"OK" if ok else "FAIL"}')
print('controls (prediction must fail):')
cfails = 0
ctrl = [(n, 2) for n in range(4, 10)] + [(6, 4), (7, 5), (8, 6), (7, 3)]
for n, r in ctrl:
    ok, full, rk = predicted(n, r)
    cfails += ok
    print(f'  N_{n}({r}): Phi^(2n-2)=I {full}, rank={rk}  {"refused OK" if not ok else "ACCEPTED (control failed)"}')
print(f'failures on N_n(n-1): {fails}; controls wrongly accepted: {cfails} of {len(ctrl)}')
