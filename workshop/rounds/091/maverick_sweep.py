"""Round 091 maverick: R_k growth of nu^k P_l for several members of each periodic indefinite key group (needs workshop/rounds/091/maverick_groups.py N first).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/091/maverick_sweep.py N K PER_GROUP [--plan]
For each group: PER_GROUP members spread through the list (first, last, evenly spaced); prints R_1..R_K and whether any k<=K has R_k == n
(a necessary sign of nu^k P ~ P[b]); same-group members must agree on the series if derived equivalent."""
import sys, os, pickle, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "089"))
from maverick_bimod import run
N, K, PER = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
g = pickle.load(open(f"/tmp/mav/g{N}.pkl", "rb"))
for gi, (key, rows) in enumerate(sorted(g.items(), key=lambda kv: len(kv[1]))):
    idx = sorted(set(round(i * (len(rows) - 1) / max(1, PER - 1)) for i in range(PER)))
    print(f"GROUP {gi} size {len(rows)} order/sign/neg {rows[0][1:]}", flush=True)
    for i in idx:
        rl = rows[i][0]
        if max(rl) > 9: print("skip", rl); continue
        s = "".join(map(str, rl)); t0 = time.time()
        if "--plan" in sys.argv: print(" ", s); continue
        A, res, _ = run(list(rl), K, verbose=False)
        R = [sum(r) for r in res]
        print(f"  {s} R_k={R[1:]} minR={min(R[1:])} n={A.n} {time.time()-t0:.0f}s", flush=True)
