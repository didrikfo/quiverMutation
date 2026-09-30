"""k(c) = n - s and hi for every core with a fit at n=13, plus every fitting centre (ambiguity). Run from repo root."""
import json, sys, os
sys.path.insert(0, "workshop/rounds/002")
from experimentalist_fit import fit
n = 13
rows = []
for l in open("workshop/rounds/002/experimentalist_census_n13.jsonl"):
    r = json.loads(l); d, s = fit(r)
    if d is None: continue
    offs = set(r["offsets"]); hi = max(offs); orbs = [set(o["held"]) for o in r["orbits"]]
    allfit = []
    for s2 in range(0, 2 * hi + 1):
        ok = sw = True, False
        ok, sw = True, False
        for ob in orbs:
            for o in ob:
                q = s2 - o
                if q in offs:
                    if q not in ob: ok = False
                    elif q != o: sw = True
        if ok and sw: allfit.append(s2)
    c = r["core"]
    rows.append((c, sum(map(int, c)), len(c), hi, n - hi, s, n - s, hi - s, allfit))
print("core digsum len hi w=n-hi s k=n-s d=hi-s allfittingcentres")
for x in sorted(rows, key=lambda x: (x[2], x[0])): print(*x)
