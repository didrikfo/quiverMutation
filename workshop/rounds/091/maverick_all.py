"""Round 091 maverick: ALL members of each periodic indefinite key group at N: min_{k=1..K} R_k, R_K, and a growth flag.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/091/maverick_all.py N K [group_index|-] [CAP] ; needs /tmp/mav/g<N>.pkl from maverick_groups.py
flag FCY? = some k <= K has R_k == n (every nu^k P_l is a single indecomposable, i.e. P_l[b]-shaped);
flag grows = R_K > R_{K/2} > R_{K/4} (linear growth envelope)."""
import sys, os, pickle, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "089"))
from maverick_bimod import run
N, K = int(sys.argv[1]), int(sys.argv[2])
g = pickle.load(open(f"/tmp/mav/g{N}.pkl", "rb"))
only = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != '-' else None
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else None   # evenly spaced sample of at most CAP members per group
for gi, (key, rows) in enumerate(sorted(g.items(), key=lambda kv: len(kv[1]))):
    if only is not None and gi != only: continue
    if CAP and len(rows) > CAP: rows = [rows[round(i * (len(rows) - 1) / (CAP - 1))] for i in range(CAP)]
    cnt = collections.Counter(); t0 = time.time(); mins = []
    for rl, o, s, neg in rows:
        if max(rl) > 9: cnt["skip"] += 1; continue
        A, res, _ = run(list(rl), K, verbose=False); R = [sum(r) for r in res]
        fcy = any(R[k] == A.n for k in range(1, K + 1))
        grows = R[K] > R[K // 2] > R[K // 4]
        cnt[("FCY" if fcy else "nonFCY", "grows" if grows else "flat")] += 1; mins.append(min(R[1:]))
    print(f"GROUP {gi} members run {len(rows)} of {len(g[key])} order/sign {rows[0][1:3]} n={N} K={K}: {dict(cnt)} min-over-members of min_k R_k = {min(mins)} {time.time()-t0:.0f}s", flush=True)
