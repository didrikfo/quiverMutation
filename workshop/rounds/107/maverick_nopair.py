"""Round 107 maverick: look for a certified-inequivalent pair at lam < 1.5, n = 8..12 (Cartan-invariant certificate, Ladkani math/0610685 Cor 3.13:
derived equivalent => C_a, C_b Z-congruent => det S, Smith form of S = C+C^T, and the discriminant form q_S agree).  Inside each lnaKeys group
(equal Phi charpoly, S signature (n-1,1), det C = 1, one real eigenvalue lam in (1,1.5)) sample members (seeded), tabulate (|det S|, SNF), and for the
first two SNF classes sharing (|det|,SNF) test q_S isomorphism against a class representative.  usage (repo root):
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_nopair.py [sample] [nmin] [nmax]"""
import sys, collections, random; sys.path.insert(0, "workshop/rounds/097")
import numpy as np, flint
from maverick_cartanZ import *
from quivermutation import coxeterTables as ct
samp = int(sys.argv[1]) if len(sys.argv) > 1 else 300
nmin = int(sys.argv[2]) if len(sys.argv) > 2 else 8; nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 12
def lam_of(n, rl):
    C = np.array(ct.lnaCartanMatrix(n, list(rl))).astype(float)
    if round(np.linalg.det(C)) != 1: return None
    ev = np.linalg.eigvalsh(C + C.T)
    if (ev < -1e-9).sum() != 1 or (abs(ev) < 1e-9).any(): return None
    w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C); big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-6]
    return big[0] if len(big) == 1 else None
random.seed(7)
for n in range(nmin, nmax + 1):
    keys = ct.lnaKeys(n); groups = collections.defaultdict(list)
    for rl, kk in keys.items(): groups[kk].append(tuple(rl))
    for kk, mem in groups.items():
        lam = lam_of(n, mem[0])
        if lam is None or lam >= 1.5: continue
        sub = random.sample(mem, min(samp, len(mem))); cls = collections.defaultdict(list)
        for rl in sub:
            S = flint.fmpz_mat([[int(x) for x in r] for r in (cartan(n, rl) + cartan(n, rl).T).tolist()])
            snf = tuple(abs(int(S.snf()[i, i])) for i in range(n)); cls[(abs(int(S.det())), snf)].append(rl)
        print("n", n, "lam %.4f" % lam, "group", len(mem), "sampled", len(sub), "(|detS|,SNF) classes", len(cls), {k[0]: len(v) for k, v in cls.items()}, flush=True)
        if len(cls) > 1: print("   CERTIFIED INEQUIVALENT by det/SNF:", [v[0] for v in list(cls.values())[:2]], flush=True)
        for k, v in cls.items():
            if k[0] <= 1 or k[0] > 3000 or len(v) < 2: continue
            reps = []
            for rl in v[:40]:
                S = cartan(n, rl) + cartan(n, rl).T; hit = False
                for r in reps:
                    try: ok, _ = disc_iso(r[1], S)
                    except RuntimeError: ok = True
                    if ok: hit = True; break
                if not hit: reps.append((rl, S))
            print("   |detS|", k[0], "q_S classes among", min(40, len(v)), ":", len(reps), [r[0] for r in reps] if len(reps) > 1 else "", flush=True)
