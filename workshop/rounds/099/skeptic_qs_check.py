"""Round 099 skeptic: independent check of q_S(90000000) ~ q_S(00000030), not using disc_iso.
Enumerates D = S^-1 Z^n / Z^n by brute force (Smith-free: columns of S^-1 mod 1 generate), then tries all bijections preserving q = x^T S x mod 2.
usage: .venv/bin/python workshop/rounds/099/skeptic_qs_check.py"""
import sys, itertools, types
sys.path.insert(0, "workshop/rounds/097"); sys.modules.setdefault("cypari", types.SimpleNamespace(pari=None))
import sympy as sp
from fractions import Fraction as F
from maverick_cartanZ import cartan
def S_of(s):
    C = cartan(10, [int(c) for c in s]); return C + C.T
def disc(S):
    Si = S.inv(); n = S.shape[0]
    gens = [tuple(F(int(Si[i, j].p), int(Si[i, j].q)) % 1 for i in range(n)) for j in range(n)]
    el = {tuple(F(0) for _ in range(n))}
    while True:
        new = set(el)
        for e in el:
            for g in gens: new.add(tuple((a + b) % 1 for a, b in zip(e, g)))
        if new == el: break
        el = new
    return sorted(el)
def q(S, x):
    n = S.shape[0]; return sum(x[i] * int(S[i, j]) * x[j] for i in range(n) for j in range(n)) % 2
def b(S, x, y):
    n = S.shape[0]; return sum(x[i] * int(S[i, j]) * y[j] for i in range(n) for j in range(n)) % 1
add = lambda x, y: tuple((a + c) % 1 for a, c in zip(x, y))
SA, SB = S_of("90000000"), S_of("00000030")
DA, DB = disc(SA), disc(SB)
print("orders", len(DA), len(DB))
print("q values A:", [str(q(SA, x)) for x in DA]); print("q values B:", [str(q(SB, x)) for x in DB])
found = 0
for perm in itertools.permutations(DB):
    m = dict(zip(DA, perm))
    if all(q(SB, m[x]) == q(SA, x) for x in DA) and all(m[add(x, y)] == add(m[x], m[y]) for x in DA for y in DA): found += 1
print("group isomorphisms preserving q:", found)
# power check: does q_S alone (no signature) separate 90000000 from 34504030 (sig (8,2))?
SC = S_of("34504030"); DC = disc(SC)
print("q values C (34504030):", [str(q(SC, x)) for x in DC], "orders", len(DC))
f2 = 0
for perm in itertools.permutations(DC):
    m = dict(zip(DA, perm))
    if all(q(SC, m[x]) == q(SA, x) for x in DA) and all(m[add(x, y)] == add(m[x], m[y]) for x in DA for y in DA): f2 += 1
print("q-preserving isos A -> C (34504030):", f2)
