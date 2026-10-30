"""Round 099 skeptic. Does q_S (discriminant form of S = C + C^T) split 90000000 vs 00000030 (n = 10, E-178/E-190, equivalent)?
HIT (refutation of q_S as a derived invariant) = q_S non-isomorphic for these two rows, or for any two rows of their key group joined
by D_10 membership (equivalent per E-190).  MISS = q_S isomorphic.  Also lists q_S classes over the whole key group (n=10).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/099/skeptic_qs_periodic.py"""
import sys, collections
sys.path.insert(0, "workshop/rounds/097")
import types; sys.modules.setdefault("cypari", types.SimpleNamespace(pari=None))
import numpy as np, sympy as sp
from quivermutation import coxeterTables as ct
from maverick_cartanZ import cartan, disc_iso, sig

N = 10
def row(s): return tuple(int(c) for c in s)
A, B = row("90000000"), row("00000030")
keys = ct.lnaKeys(N)
print("keys equal:", keys[A] == keys[B])
CA, CB = cartan(N, list(A)), cartan(N, list(B))
SA, SB = CA + CA.T, CB + CB.T
print("det", int(SA.det()), int(SB.det()), "sig", sig(SA), sig(SB))
print("q_S iso A~B:", disc_iso(SA, SB))
# whole key group: partition by q_S iso class (representative comparison)
grp = [rl for rl, k in keys.items() if k == keys[A]]
print("key group size", len(grp))
reps = []; cnt = collections.Counter(); member = {}
for rl in grp:
    C = cartan(N, list(rl)); S = C + C.T
    for i, (S0, r0) in enumerate(reps):
        if int(S.det()) == int(S0.det()) and sig(S) == sig(S0) and disc_iso(S0, S)[0]:
            cnt[i] += 1; member[rl] = i; break
    else:
        reps.append((S, rl)); cnt[len(reps)-1] += 1; member[rl] = len(reps)-1
print("q_S classes in key group:", [(i, reps[i][1], int(reps[i][0].det()), sig(reps[i][0]), cnt[i]) for i in range(len(reps))])
print("class of A:", member[A], "class of B:", member[B])
