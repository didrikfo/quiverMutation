"""Round 097 maverick: controls for maverick_cartanZ.disc_iso: (1) S vs P^T S P for random unimodular P (must be isomorphic),
(2) S vs S with one pair of vertices permuted, (3) all members of a group: number of distinct discriminant-form classes.
usage (repo root): .venv/bin/python workshop/rounds/097/maverick_control.py"""
import sys, random, sympy as sp
sys.path.insert(0, "workshop/rounds/097")
from maverick_cartanZ import *
N = 11
a = (2, 2, 2, 2, 2, 3, 0, 3, 0); b = (2, 2, 4, 0, 0, 2, 2, 3, 0)
for rl in (a, b):
    C = cartan(N, rl); S = C + C.T
    rng = random.Random(1); P = sp.eye(N)
    for _ in range(40):
        i, j = rng.sample(range(N), 2); E = sp.eye(N); E[i, j] = rng.choice([-1, 1, 2]); P = P * E
    print(rl, "random congruence:", disc_iso(S, P.T * S * P)[0], "det P", P.det())
    perm = list(range(N)); perm[0], perm[1] = perm[1], perm[0]
    print(" perm:", disc_iso(S, S.extract(perm, perm))[0])
Sa = cartan(N, a); Sa = Sa + Sa.T; Sb = cartan(N, b); Sb = Sb + Sb.T
print("a vs b", disc_iso(Sa, Sb)[0], "b vs a", disc_iso(Sb, Sa)[0])
# invariant: p-rank / group structure
for S in (Sa, Sb):
    d, el, gens, bil, q, add = disc_form(S)
    from collections import Counter
    print("order", d, "gens", len(gens), "q-value histogram", sorted(Counter((order(x, add), (sum(1 for _ in [0]) and 0)) for x in el).items()))
print("--- independent look: multiset of (order, x^T S x mod 2) over all elements")
from collections import Counter
for S in (Sa, Sb):
    d, el, gens, bil, q, add = disc_form(S); n = N
    Sf = [[Fraction(int(S[i, j])) for j in range(n)] for i in range(n)]
    qq = lambda x: sum(x[i]*Sf[i][j]*x[j] for i in range(n) for j in range(n)) % 2
    print(sorted(Counter((order(x, add), qq(x)) for x in el).items()))
