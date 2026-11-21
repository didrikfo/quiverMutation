# nu^m on every uniserial module [s,e] (top s, socle e) of a linear Nakayama algebra; report returns M[c] and their shifts.
# Adapted from the round-134 referee's scratch script unis.py (maverick). Usage: scholar_uniserial.py REL m1,m2 [--prime p]
# A return is counted only when H^* has the dimension vector of [s,e] in a single degree; since nu is an equivalence the
# object is indecomposable, hence that cohomology module is a uniserial with top s and length e-s+1, i.e. [s,e].
import sys, time
R = "/home/user/quiverMutation/workshop/rounds"
sys.path.insert(0, R + "/134"); sys.path.insert(0, R + "/089"); sys.path.insert(0, R + "/087")
if "--prime" in sys.argv:
    i = sys.argv.index("--prime"); p = int(sys.argv[i + 1]); del sys.argv[i:i + 2]
    import scholar_rankgrowth, maverick_bimod, scholar_simple_nu
    scholar_rankgrowth.P = maverick_bimod.P = scholar_simple_nu.P = p
from maverick_bimod import resolve, nu_bimod
from scholar_rankgrowth import Alg, Cx
from scholar_simple_nu import cohomology
def res(A, s, e):
    ivs = [(s, e)]
    while True:
        a, b = ivs[-1]; nb = a + A.p[a] - 1
        if b + 1 > nb: break
        ivs.append((b + 1, nb))
    X = Cx()
    for k, (a, _) in enumerate(ivs): X.add(k, -k, a)
    for k in range(1, len(ivs)): X.setd(k, k - 1, 1)
    return X
print("P =", __import__("scholar_rankgrowth").P, flush=True)
rel = sys.argv[1]; ms = [int(x) for x in sys.argv[2].split(",")]; A = Alg([int(c) for c in rel]); n = A.n; L = resolve(A); t0 = time.time()
out = {m: [] for m in ms}
for s in range(1, n + 1):
    for e in range(s, s + A.p[s]):
        X = res(A, s, e); H0 = cohomology(A, X); assert H0 == {(i, 0): 1 for i in range(s, e + 1)}, (s, e, H0)
        k = 0
        for m in ms:
            while k < m: X = nu_bimod(A, L, X)[0]; k += 1
            H = cohomology(A, X); ts = {t for (_, t) in H}
            if len(ts) == 1 and H == {(i, list(ts)[0]): 1 for i in range(s, e + 1)}: out[m].append((s, e, -list(ts)[0]))
for m in ms:
    sh = sorted({c for *_, c in out[m]})
    print(rel, "m", m, "returns", len(out[m]), "shifts", sh, {c: [f"[{s},{e}]" for s, e, cc in out[m] if cc == c] for c in sh}, f"{time.time()-t0:.0f}s", flush=True)
