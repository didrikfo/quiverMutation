"""Round 131 scholar: objects-only rank growth of nu^k P_l (E-246's minimal-bimodule method, reused unchanged from
workshop/rounds/089/maverick_bimod.py) on the radical-power LNAs that arXiv:2112.15587v2 places in Fuchsian type.
Prints min_k R_k, every k with R_k = n (R_k = n is necessary for S^k = [b], E-245), and R_k at multiples of the
Phi-order.  R_k > n at every k <= K means S^k is not a shift on projectives for k <= K.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/131/scholar_rank.py RELSTRING K ORDER"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "089"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "087"))
from maverick_bimod import run

rel, K, order = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
A, res, cur = run([int(c) for c in rel], K, verbose=False)
n = A.n
R = [sum(r) for r in res]
print(f"{rel} n={n} K={K} Phi-order={order} time={time.time()-t0:.0f}s")
print("min R_k (k>=1):", min(R[1:]), "at k =", [k for k in range(1, K + 1) if R[k] == min(R[1:])][:8])
print("k with R_k = n:", [k for k in range(1, K + 1) if R[k] == n])
print("R_k at multiples of the order:", {m: R[m] for m in range(order, K + 1, order)})
print("R_k series:", R)
