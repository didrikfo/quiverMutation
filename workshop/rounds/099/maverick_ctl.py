"""Round 099 maverick: control for maverick_iso140.test: C_b = P C_a P^T for a random P in GL_12(Z) (small entries) must be found congruent.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/099/maverick_ctl.py"""
import sys; sys.path.insert(0, "workshop/rounds/099")
import maverick_iso140 as M
import numpy as np, random, sympy as sp
a = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0)
Ca = M.cartan(12, a); rng = random.Random(1)
P = sp.eye(12)
for _ in range(25):
    i, j = rng.sample(range(12), 2); E = sp.eye(12); E[i, j] = rng.choice([-1, 1]); P = E * P
Cb = P * Ca * P.T
M.cartan0 = M.cartan
M.cartan = lambda N, rl: Cb if rl == "B" else M.cartan0(N, rl)
for slack in (0.5, 30):
    r = M.test(a, "B", slack); print(slack, {k: v for k, v in r.items() if k != "P"})
