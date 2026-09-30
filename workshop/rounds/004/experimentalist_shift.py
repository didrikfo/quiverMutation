"""Centre s and shortfall d of the 12 E-060 cores by n (13 census, 14 sample, 15, 16 files).
  .venv/bin/python workshop/rounds/004/experimentalist_shift.py FILE.jsonl ...
Prints core n s d k=n-s pairs-flag; fit rule of round 002 (d = hi - s, |d| <= 6)."""
import json, sys
sys.path.insert(0, "workshop/rounds/002")
from experimentalist_fit import fit
for p in sys.argv[1:]:
    for l in open(p):
        r = json.loads(l); res = fit(r)
        if res[0] is None or res[1] is None:
            print(r["core"], r["n"], "nofit"); continue
        d, s = res
        print(r["core"], r["n"], "s", s, "hi", max(r["offsets"]), "k", r["n"] - s, "d", max(r["offsets"]) - s,
              "orbits", len(r["orbits"]), "allclosed", all(o["closed"] for o in r["orbits"]))
