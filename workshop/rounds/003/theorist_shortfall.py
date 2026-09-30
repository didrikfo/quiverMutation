"""Test: is the reflection shortfall d(c) related to the H-020 head/tail?

  .venv/bin/python workshop/rounds/003/theorist_shortfall.py
Inputs: round-002 census (orbits), theorist_slides_n13.jsonl (slides).
d, s from experimentalist_fit.fit; hi = last offset; head/tail from batch._headAndTail.
"""
import json, os, sys, collections
sys.path.insert(0, os.getcwd())
sys.path.insert(0, "workshop/rounds/002")
import batch
from experimentalist_fit import fit
cen = {}
for l in open("workshop/rounds/002/experimentalist_census_n13.jsonl"):
    r = json.loads(l); cen[r["core"]] = r
sl = {}
for l in open("workshop/rounds/003/theorist_slides_n13.jsonl"):
    r = json.loads(l); sl[r["core"]] = r["slide"]
rows = []
for c, r in cen.items():
    d, s = fit(r)
    if d is None: continue
    slide = sl[c]; h, t, interior = batch._headAndTail(slide)
    hi = max(r["offsets"]); 
    rows.append((c, slide, h, t, interior, d, s, hi, min(r["offsets"])))
print("core slide head tail interior d s hi lo  |t-h| s-hi")
cnt = collections.Counter()
for c, slide, h, t, i, d, s, hi, lo in rows:
    print(c, slide, h, t, i, d, s, hi, lo, abs(t-h), hi-s)
    cnt[(d == abs(t-h), d == max(t-h,0), d == max(h-t,0))] += 1
print(cnt)

print("\n== test s == hi + head - tail (signed shortfall d = hi - s = tail - head)")
def cls(slide): return "allI" if set(slide)=={"i"} else ("allO" if set(slide)=={"o"} else "mixed")
res = collections.defaultdict(list)
for c, slide, h, t, i, d, s, hi, lo in rows:
    res[cls(slide)].append((c, s == hi + h - t))
for k, v in res.items():
    print(k, sum(x for _, x in v), "/", len(v), "fail:", [c for c, x in v if not x])
print("cores without a reflection fit, slides:")
for c in cen:
    if fit(cen[c])[0] is None: print(" ", c, sl[c], [o["held"] for o in cen[c]["orbits"]])
