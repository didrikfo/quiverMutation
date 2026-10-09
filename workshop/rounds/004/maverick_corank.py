"""Signature (pos, neg, zero) of the symmetrised Euler form C+C^T of every LNA of length N, by quipu status.
usage: python workshop/rounds/004/maverick_corank.py N"""
import sys, collections, numpy as np
from quivermutation import coxeterTables as ct, nakayama as nk, invariants as inv
n = int(sys.argv[1]); st = ct.lnaStatus(n); tab = collections.Counter()
for rl, s in st.items():
    C = np.array(inv.cartanMatrix(nk.LinearNakayamaAlgebra(n, list(rl))), dtype=float)
    ev = np.linalg.eigvalsh(C + C.T)
    tab[(s, int((ev > 1e-8).sum()), int((ev < -1e-8).sum()), int((abs(ev) <= 1e-8).sum()))] += 1
for k, v in sorted(tab.items()): print(k, v)
