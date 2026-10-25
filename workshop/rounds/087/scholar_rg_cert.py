"""Round 087 scholar: object-level certificate for nu^a P_l ~ P_sigma(l)[b_l].  Runs nu a times on every P_l and prints, for each l, the
minimal complex of nu^a P_l as a list of (degree, vertex); a single summand for every l, with the same degree for all, is the
projective-level form of S^a = [b] (it is necessary for S^a = [b] as functors, not sufficient).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/087/scholar_rg_cert.py RELSTRING a"""
import sys
sys.path.insert(0, "workshop/rounds/087")
import scholar_rankgrowth as rg
rl = [int(c) for c in sys.argv[1]]; a = int(sys.argv[2])
A, res, cur = rg.run(rl, a, 400000, verbose=False)
print("totals", [sum(r) for r in res])
for l, X in cur.items(): print(f"nu^{a} P_{l} =", sorted(X.vert.values()))
