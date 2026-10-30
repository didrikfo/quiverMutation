"""Round 099 maverick: the opposite-algebra convention: is C_b^T Z-congruent to C_a (P C_a P^T = C_b^T)?  for the open 140 pair and for a within-quipu pair (control).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/099/maverick_transp.py"""
import sys; sys.path.insert(0, "workshop/rounds/099")
import maverick_iso140 as M
a = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); b = (2, 2, 3, 0, 3, 0, 2, 2, 3, 0); a2 = (2, 2, 2, 3, 0, 0, 3, 0, 3, 0)
M.cartan0 = M.cartan
M.cartan = lambda N, rl: M.cartan0(N, rl[1]).T if rl[0] == "T" else M.cartan0(N, rl)
for name, x, y in (("open pair a vs b^T", a, ("T", b)), ("control a vs a2 (within quipu) ^T", a, ("T", a2)), ("a vs a^T", a, ("T", a))):
    r = M.test(x, y, 0.5); print(name, {k: v for k, v in r.items() if k != "P"})
