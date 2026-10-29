"""Round 097 maverick: Phi-equivariant refinement for the n = 12 pair that q_S does not split: for each irreducible factor f of the Coxeter polynomial
(and each prime-power-free product), the primitive sublattice K_f = ker f(Phi) in Z^n with Gram S|K_f; compare rank, signature, det and discriminant form
(disc_iso).  A difference in any piece is a Cartan-level NO for C-congruence (P maps K_f(a) onto K_f(b) isometrically).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/097/maverick_pieces.py"""
import sys; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
N = 12; x = sp.symbols('x')
a = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); b = (2, 2, 3, 0, 3, 0, 2, 2, 3, 0)
def pieces(rl):
    C = cartan(N, rl); S = C + C.T; Ph = phi(C); out = {}
    cp = sp.Poly(Ph.charpoly(x).as_expr(), x); fs = [g for g, e in sp.factor_list(cp.as_expr())[1]]
    for f in fs + [sp.expand((x + 1) ** 2)]:
        fm = sp.zeros(N, N); Pw = sp.eye(N)
        for c in reversed(sp.Poly(f, x).all_coeffs()): fm += c * Pw; Pw = Pw * Ph
        Km = pari.matkerint(pari.matrix(N, N, [int(v) for v in fm])); k = int(pari.matsize(Km)[1])
        Bm = sp.Matrix(N, k, [int(Km[i, j]) for i in range(N) for j in range(k)])
        out[str(f)] = (Bm.T * S * Bm)
    return out
pa, pb = pieces(a), pieces(b)
for f in pa:
    Ga, Gb = pa[f], pb[f]; k = Ga.shape[0]
    print("piece", f, "rank", k, "det", Ga.det(), Gb.det(), "sig", sig(Ga), sig(Gb), end=" ")
    if Ga.det() != 0 and Ga.det() == Gb.det() and sig(Ga) == sig(Gb):
        print("q iso:", disc_iso(Ga, Gb)[0])
    else: print("(degenerate or det/sig differ)")
