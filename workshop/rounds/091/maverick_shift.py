"""Round 091 maverick: for FCY key groups, at the first k <= K with R_k = n check that nu^k P_l = P_sigma(l)[b_l]:
print sigma (is it id?) and whether the shift b_l is the same for all l.  Evenly spaced sample of M members per group.
usage: timeout 10m .venv/bin/python workshop/rounds/091/maverick_shift.py N K group_index M ; needs /tmp/mav/g<N>.pkl"""
import sys, os, pickle, collections
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "089"))
import maverick_bimod as mb
from scholar_rankgrowth import Alg
N, K, gi, M = map(int, sys.argv[1:5])
g = pickle.load(open(f"/tmp/mav/g{N}.pkl", "rb"))
key, rows = sorted(g.items(), key=lambda kv: len(kv[1]))[gi]
rows = [rows[round(i * (len(rows) - 1) / (M - 1))] for i in range(M)] if len(rows) > M else rows
summ = collections.Counter()
for rl, o, s, neg in rows:
    A = Alg(list(rl)); levels = mb.resolve(A)
    cur = {}
    for l in range(1, A.n + 1):
        X = mb.Cx(); X.vert[0] = (0, l); cur[l] = X
    first = None
    for k in range(1, K + 1):
        cur = {l: mb.nu_bimod(A, levels, X)[0] for l, X in cur.items()}
        if all(len(X.vert) == 1 for X in cur.values()): first = k; break
    if first is None: print("member", "".join(map(str, rl)), "no k <= K with R_k = n"); summ["none"] += 1; continue
    img = {l: list(X.vert.values())[0] for l, X in cur.items()}   # (degree, projective index)
    sigma = {l: v[1] for l, v in img.items()}; shifts = {v[0] for v in img.values()}
    ident = all(sigma[l] == l for l in sigma); perm = sorted(sigma.values()) == list(range(1, A.n + 1))
    print("member", "".join(map(str, rl)), f"first k={first} sigma_is_id={ident} sigma_is_perm={perm} shifts={sorted(shifts)} sigma={[sigma[l] for l in sorted(sigma)]}", flush=True)
    summ[(first, ident, perm, len(shifts) == 1)] += 1
print("SUMMARY n=%d group %d (size %d, order/sign %s): (first k, sigma=id, sigma perm, common shift) -> count: %s" % (N, gi, len(g[key]), rows[0][1:3], dict(summ)))
