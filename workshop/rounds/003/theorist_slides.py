"""Slides (i/o/? per offset) at n = 13 for the 139 cores of the round-002 census.

Verdict of each recorded orbit is taken at its least offset (an orbit has one
verdict, E-052), via batch._verdictFor under the reduced walk.  Resumable.

  .venv/bin/python workshop/rounds/003/theorist_slides.py K/T OUT.jsonl [CENSUS.jsonl [N]]
"""
import json, os, sys
sys.path.insert(0, os.getcwd())
import batch
from quivermutation import freeMoves as fm

k, t = (int(x) for x in sys.argv[1].split("/"))
out = sys.argv[2]
censusPath = sys.argv[3] if len(sys.argv) > 3 else "workshop/rounds/002/experimentalist_census_n13.jsonl"
N = int(sys.argv[4]) if len(sys.argv) > 4 else 13
census = [json.loads(l) for l in open(censusPath)][k::t]
done = set()
if os.path.exists(out):
    done = {json.loads(l)["core"] for l in open(out)}
with open(out, "a") as h:
    for rec in census:
        if rec["core"] in done:
            continue
        v = {}
        for orb in rec["orbits"]:
            o = min(orb["held"])
            r = batch._verdictFor(N, rec["core"], o, orbitLimit=1500000,
                                  joinLimit=1500000, free=fm.REDUCED)["verdict"]
            for x in orb["held"]:
                v[x] = r
        L = {"inside": "i", "outside": "o", "undecided": "?"}
        slide = "".join(L[v[x]] if x in v else "." for x in range(max(v) + 1))
        h.write(json.dumps(dict(core=rec["core"], slide=slide)) + "\n"); h.flush()
