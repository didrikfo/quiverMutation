"""Round 074 response: (1) direct Cartan-congruence check of the P returned by isometry (P C_a P^T = C_b, i.e. Q = P^T gives Q^T C_a Q = C_b);
(2) the artificial non-Cartan NO control (G, M) vs (G, -M); (3) search for a real NO between two definite Cartan data of equal key (Coxeter
key from ct.lnaKeys) and equal detG at n = 5..9.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/074/maverick_controls.py [part: 1|2|3] """
import sys, time, collections
sys.path.insert(0, "workshop/rounds/074")
import numpy as np
from maverick_isometry import data, isometry
from quivermutation import coxeterTables as ct, nakayama as nk
part = sys.argv[1] if len(sys.argv) > 1 else "all"
t0 = time.time()
Cm = lambda rl, N: np.array(ct.lnaCartanMatrix(N, tuple(rl)), dtype=np.int64)
def cartanOK(P, a, b, N):
    P = np.array(P.tolist(), dtype=np.int64)
    return bool((P @ Cm(a, N) @ P.T == Cm(b, N)).all()), bool((P.T @ Cm(b, N) @ P == Cm(a, N)).all())
rd = lambda s: tuple(int(c) for c in s)
if part in ("all", "1"):
    N = 10; keys = ct.lnaKeys(N); tot = collections.Counter(); groups = collections.defaultdict(list)
    for rl in nk.allRelationLengths(N):
        C = Cm(rl, N).astype(float); Ci = np.linalg.inv(C)
        if np.linalg.eigvalsh(Ci + Ci.T).min() < 1e-9: continue
        Phi = np.rint(-Ci.T @ C).astype(np.int64); I = np.eye(N, dtype=np.int64); M = Phi.copy(); o = 0
        for m in range(1, 201):
            if (M == I).all(): o = m; break
            M = M @ Phi
        if o: groups[keys[rl]].append(rl)
    for k, mem in groups.items():
        a = mem[0]; Ga, Ma, _ = data(a, N)
        for b in mem[1:]:
            Gb, Mb, _ = data(b, N); r, P, _ = isometry(Ga, Ma, Gb, Mb)
            if r == "YES":
                f, bw = cartanOK(P, a, b, N); tot[("YES", "P C_a P^T=C_b", f, "P^T C_b P=C_a", bw)] += 1
            else: tot[r] += 1
    print("part1 (n=10 survey, member b vs rep a):", dict(tot), f"{time.time()-t0:.0f}s", flush=True)
if part in ("all", "2"):
    G, M, _ = data(rd("90000000"), 10); r, _, nodes = isometry(G, M, G, -M)
    print("part2 (G,M) vs (G,-M), 90000000 n=10 [non-Cartan artificial]:", r, "nodes", nodes, flush=True)
if part in ("all", "3"):
    for N in range(4, 10):
        keys = ct.lnaKeys(N); grp = collections.defaultdict(list)
        for rl in nk.allRelationLengths(N):
            C = Cm(rl, N).astype(float); Ci = np.linalg.inv(C)
            if np.linalg.eigvalsh(Ci + Ci.T).min() < 1e-9: continue
            grp[(keys[rl], round(np.linalg.det(Ci + Ci.T) * round(np.linalg.det(C)) ** 2))].append(rl)
        c = collections.Counter(); nos = []
        for k, mem in grp.items():
            a = mem[0]; Ga, Ma, _ = data(a, N)
            for b in mem[1:]:
                Gb, Mb, _ = data(b, N); r, P, _ = isometry(Ga, Ma, Gb, Mb); c[r] += 1
                if r == "NO": nos.append(("".join(map(str, a)), "".join(map(str, b))))
                if r == "YES" and not cartanOK(P, a, b, N)[0]: c["YES-but-not-Cartan"] += 1
            if time.time() - t0 > 500: break
        print(f"part3 n={N}: key+detG groups {len(grp)} (members vs first) {dict(c)} NO pairs {nos[:8]}", f"{time.time()-t0:.0f}s", flush=True)
