"""Round 134 scholar: two different shifts under nu^m exclude fractional Calabi-Yau (H-025, E-339 blocker).
For each simple S_l and each projective P_l, test whether nu^m X ~ X[c] (exact: for S_l, cohomology is S_l in one
degree; for P_l, the minimal complex is P_l in one degree).  If two objects return with shifts c1 != c2 at the same m,
then nu^a ~ [b] (as functors, a >= 1) is impossible: it would give X[a c] ~ nu^{m a} X ~ X[m b] for both X,
so a c1 = m b = a c2.  On a hypersurface-type control all returns must have one shift (S^h = [h + 2]).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/134/scholar_two_slopes.py RELSTRING m"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "089")); sys.path.insert(0, os.path.join(HERE, "..", "087"))
from maverick_bimod import resolve, nu_bimod
from scholar_rankgrowth import Alg, Cx
from scholar_simple_nu import simple_resolution, cohomology

rel, m = sys.argv[1], int(sys.argv[2]); A = Alg([int(c) for c in rel]); n = A.n; t0 = time.time()
L = resolve(A); returns = []
for l in range(1, n + 1):
    for kind in ("S", "P"):
        if kind == "S": X = simple_resolution(A, l)
        else: X = Cx(); X.vert[0] = (0, l)
        for _ in range(m): X = nu_bimod(A, L, X)[0]
        if kind == "S":
            H = cohomology(A, X)
            if len(H) == 1 and list(H.values()) == [1] and list(H)[0][0] == l: returns.append((kind, l, -list(H)[0][1]))
        elif len(X.vert) == 1 and list(X.vert.values())[0][1] == l:
            returns.append((kind, l, -list(X.vert.values())[0][0]))
shifts = sorted(set(c for _, _, c in returns))
print(f"{rel} n={n} m={m}: {len(returns)} of {2*n} simples/projectives return as X[c]; shifts {shifts}  ({time.time()-t0:.0f}s)")
for c in shifts: print(f"  c = {c}: " + " ".join(f"{k}_{l}" for k, l, cc in returns if cc == c))
print("VERDICT:", "two shifts at one m => not fractionally CY (no a >= 1, b with nu^a ~ [b])" if len(shifts) > 1
      else "single shift: consistent with FCY (no conclusion)" if shifts else "no returns: no conclusion")
