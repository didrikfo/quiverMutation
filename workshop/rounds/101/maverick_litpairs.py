"""Round 101 maverick: literature-certified derived-equivalent radical-power LNA pairs N(n,l+1)~N(n,l) (Ueda 2302.02880 Cor 1.3, Cor 1.6; Ladkani 0911.5137 Cor 1.2)
at n = 9..12, as positive controls for the Cartan-level NO test (q_S, signature, det, char poly of Phi, Z[Phi]-module search, rational Phi-conjugacy).
N(n,l) = KA_n/rad^l has digit string "l"*(n-l)+"0"*(l-2) (length n-2).  Every pair must pass every test; a NO would void the test (or the literature).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/101/maverick_litpairs.py"""
import sys; sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
def rl(n, l): return tuple([l] * (n - l) + [0] * (l - 2))
pairs = []  # (n, l1, l2, source)
for q in (1, 2):
    for r2 in range(0, 12):   # r = r2/2, p = 2
        n = 6 * q + r2; l = 3 * q + r2   # n=p(p+1)q+p(p-1)r = 6q + 2r, l = (p+1)q+pr = 3q + 2r (r = r2/2)
        if 9 <= n <= 12: pairs.append((n, l + 1, l, "Ueda Cor1.3 p=2 q=%d r=%s" % (q, r2 / 2)))
pairs.append((12, 5, 4, "Ueda Cor1.3 p=3 q=1 r=0"))
for (p, q) in ((2, 4), (4, 2), (2, 5), (5, 2), (3, 3)):
    n = p * q + 1
    if 9 <= n <= 12 and p <= q: pairs.append((n, q + 1, p + 1, "Ueda Cor1.6 p=%d q=%d" % (p, q)))
for m, k in ((2, 5), (2, 6), (3, 4), (2, 4.5)):
    if k == int(k) and 9 <= m * k <= 12: pairs.append((int(m * k), m + 1, int(k) + 1, "Ladkani Cor1.2 m=%d n=%d" % (m, k)))
seen = set(); out = []
for n, l1, l2, src in pairs:
    key = (n, max(l1, l2), min(l1, l2))
    if key in seen or l1 == l2 or l1 >= n or l2 >= n: continue
    seen.add(key); out.append((n, l1, l2, src))
for n, l1, l2, src in sorted(out):
    a, b = rl(n, l1), rl(n, l2); Ca, Cb = cartan(n, a), cartan(n, b); Sa, Sb = Ca + Ca.T, Cb + Cb.T
    Pa, Pb = phi(Ca), phi(Cb); x = sp.symbols("x")
    cp = lambda P: sp.factor_list(P.charpoly(x).as_expr())
    res = dict(n=n, pair="N(%d,%d)~N(%d,%d)" % (n, l1, n, l2), src=src, sigA=sig(Sa), sigB=sig(Sb), detA=int(Sa.det()), detB=int(Sb.det()),
               charpoly_eq=(sp.expand(Pa.charpoly(x).as_expr() - Pb.charpoly(x).as_expr()) == 0), PhiA_PhiB_modules=None, qS=None, Zcong=None)
    try: res["qS"] = disc_iso(Sa, Sb)[0] if abs(res["detA"]) > 1 else "vacuous(|det|<=1)"
    except Exception as e: res["qS"] = "err:" + str(e)[:40]
    try:
        r = module_iso_search(Pa, Pb, 3000)
        res["PhiA_PhiB_modules"] = r[:3]
    except Exception as e: res["PhiA_PhiB_modules"] = "err:" + str(e)[:40]
    print(res, flush=True)
