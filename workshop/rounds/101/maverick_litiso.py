"""Round 101 maverick: run the E-269 exhaustive Phi-commutant Z-congruence test (workshop/rounds/099/maverick_iso140.py, test()) on literature-certified equivalent
radical-power pairs with S of signature (n-1,1) and det S != 0 (the only ones the method applies to). A NO here would void the test.  Also a negative-control pair
(same n, different class) is NOT available with equal det/sig, so power here is YES-only.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/101/maverick_litiso.py"""
import sys; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/101")
import maverick_iso140 as M
import numpy as np
def pq(S, P):   # as in 099 but threshold 1.05 (Salem number of these pairs is 1.18..1.26, not > 1.5)
    N = M.N; w, V = np.linalg.eig(P)
    i = [k for k in range(N) if abs(w[k].imag) < 1e-9 and w[k].real > 1.05][0]; j = [k for k in range(N) if abs(w[k].imag) < 1e-9 and abs(w[k].real - 1 / w[i].real) < 1e-9][0]
    p, q = V[:, i].real, V[:, j].real; c = p @ S @ q
    return p, q / c, w[i].real
M.pq = pq
def rl(n, l): return tuple([l] * (n - l) + [0] * (l - 2))
for n, l1, l2 in ((10, 8, 7), (11, 9, 8), (12, 10, 9)):
    M.N = n
    for slack in (0.5, 30):
        try:
            r = M.test(rl(n, l1), rl(n, l2), slack)
            print("N(%d,%d)~N(%d,%d) slack %s" % (n, l1, n, l2, slack), {k: v for k, v in r.items() if k != "P"}, "P found" if r["P"] else "", flush=True)
        except Exception as e: print("N(%d,%d)~N(%d,%d)" % (n, l1, n, l2), "ERR", repr(e)[:200], flush=True)
