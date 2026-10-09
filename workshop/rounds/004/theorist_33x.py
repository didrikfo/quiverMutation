"""33x cores across n: hi, fitted centre s, k = n-s, d = hi-s.  Run from repo root after the census files exist:
  for n in 14 15 16 17; do .venv/bin/python workshop/rounds/002/experimentalist_census.py $n --cores 33,333,334,335,336 --out workshop/rounds/004/theorist_33x_n$n.jsonl; done   (~1 min each, parallel ok)
"""
import json, sys
sys.path.insert(0, "workshop/rounds/002")
from experimentalist_fit import fit
for n in (13, 14, 15, 16, 17):
    f = "workshop/rounds/004/theorist_33x_n%d.jsonl" % n if n > 13 else "workshop/rounds/002/experimentalist_census_n13.jsonl"
    for l in open(f):
        r = json.loads(l)
        if r["core"] in ("33", "333", "334", "335", "336"):
            d, s = fit(r); hi = max(r["offsets"])
            print(n, r["core"], "hi", hi, "s", s, "k", None if s is None else n - s, "d", None if s is None else hi - s,
                  [o["held"] for o in r["orbits"]] if n > 13 else "")
