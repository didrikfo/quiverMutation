"""Round 110 skeptic: blocking power check for H-024. Pairs inequivalent BY CONSTRUCTION (same key => same Phi charpoly, same |det S|, signature (n-1,1); different q_S)
inside E-269's domain; run E-269's test (workshop/rounds/099/maverick_iso140.py test()) on them. Expected: congruent False (NO). Controls: within-quipu pairs (YES) and self pairs.
Independent q_S certificate: histogram of x^T S x mod 2 over the discriminant group (own code, not disc_iso) AND maverick's disc_iso.
Round 110 variant: reports nvec/nsol at slack 0.5 AND at slack S (default 5) for the cross pair and the within-quipu control.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/110/skeptic_qsno.py N [max_groups] [slack]"""
import numpy as np, sys, collections, itertools, time
sys.path.insert(0, "workshop/rounds/097"); sys.path.insert(0, "workshop/rounds/099")
import maverick_iso140 as iso
from maverick_cartanZ import *
from fractions import Fraction
N = int(sys.argv[1]); maxg = int(sys.argv[2]) if len(sys.argv) > 2 else 4; slack = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
iso.N = N
def hist(S):
    d, E, G, bil, q, add = disc_form(S); n = S.shape[0]
    Sf = [[Fraction(int(S[i, j])) for j in range(n)] for i in range(n)]
    h = collections.Counter()
    for x in E: h[(sum(x[i] * Sf[i][j] * x[j] for i in range(n) for j in range(n)) % 2)] += 1
    return d, tuple(sorted(h.items()))
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
done = 0
for k, v in sorted(groups.items(), key=lambda kv: len(kv[1])):
    q = collections.defaultdict(list)
    for rl in v:
        r = quipuForAlmostSeparateLNA(N, list(rl))
        if r is not None: q[r].append(tuple(rl))
    if len(q) < 2: continue
    qs = list(q); a = q[qs[0]][0]; b = q[qs[1]][0]; a2 = q[qs[0]][1] if len(q[qs[0]]) > 1 else None
    Ca, Cb = cartan(N, a), cartan(N, b); Sa, Sb = Ca + Ca.T, Cb + Cb.T
    if sig(Sa) != (N - 1, 1) or sig(Sb) != (N - 1, 1) or int(Sa.det()) != int(Sb.det()) or phi(Ca).charpoly() != phi(Cb).charpoly(): continue
    ha, hb = hist(Sa), hist(Sb); io = disc_iso(Sa, Sb)[0]
    t = time.time(); lam = max(abs(np.linalg.eigvals(np.array(phi(Ca).tolist(), dtype=float))))
    out = {}
    for sl in (0.5, slack):
        try:
            r = iso.test(a, b, sl); c = iso.test(a, a2, sl) if a2 else None
            out[sl] = ((r["congruent"], r["nvec"], r["nsol"]), (c["congruent"], c["nvec"], c["nsol"]) if c else "n/a")
        except AssertionError as e:
            out[sl] = ("assert", str(e)[:60])
    print(f"N={N} size {len(v)} lam={lam:.4f} hist_differs={ha != hb} cross(cong,nvec,nsol)@0.5={out[0.5][0]} ctrl@0.5={out[0.5][1]} cross@{slack}={out[slack][0]} ctrl@{slack}={out[slack][1]} {time.time()-t:.0f}s", flush=True)
    done += 1
    if done >= maxg: break
