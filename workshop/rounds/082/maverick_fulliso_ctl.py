"""Round 082: controls for maverick_fulliso.py. (1) the 9 reps' Cartan matrices are pairwise distinct and the P found is not the identity;
(2) NO control: (G,M) vs (G,M^5) at n=10 (50505000) and n=11 (rep 0): P M = M^5 P would move the negative eigenplane (angle 1/14 -> 5/14); the search,
run with each side's OWN majorant, must return NO (exhausted, nodes counted).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/082/maverick_fulliso_ctl.py"""
import sys, time
sys.path.insert(0, "workshop/rounds/074"); sys.path.insert(0, "workshop/rounds/082")
import numpy as np
from maverick_isometry import data
from maverick_fulliso import majorant, isometryH
from maverick_eigsig import isoClasses15
from quivermutation import coxeterTables as ct
t0 = time.time(); rd = lambda s: tuple(int(c) for c in s)
def twist(a, N, k, negb):
    Ga, Ma, _ = data(a, N); Mb = np.linalg.matrix_power(Ma, k)
    Ha = majorant(Ga, Ma, {round(negb[0], 4), round(1 - negb[0], 4)} if False else NEG[N]); Hb = majorant(Ga, Mb, negb)
    return isometryH(Ga, Ma, Ga, Mb, Ha, Hb)[:3]
NEG = {10: {0.0556, 0.9444}, 11: {0.0714, 0.9286}}
print("n10 50505000: (G,M) vs (G,M^5)  neg angle 5/18=0.2778 ->", twist(rd("50505000"), 10, 5, {0.2778, 0.7222}))
mem, reps, find = isoClasses15()
Cs = [ct.lnaCartanMatrix(11, r) for r in reps]
print("distinct Cartan matrices among 9 reps:", len({tuple(map(tuple, c)) for c in Cs}))
print("n11 rep0: (G,M) vs (G,M^5), neg angle 5/14 ->", twist(reps[0], 11, 5, {0.3571, 0.6429}), f"{time.time()-t0:.0f}s")
