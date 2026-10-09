"""Rule test: fit s == hi + head - tail on mixed slides.  args: census slides
  .venv/bin/python workshop/rounds/003/theorist_rule.py workshop/rounds/003/theorist_census_n14_sample.jsonl workshop/rounds/003/theorist_slides_n14_sample.jsonl
"""
import json, os, sys
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "workshop/rounds/002")
import batch
from experimentalist_fit import fit
cen = [json.loads(l) for l in open(sys.argv[1])]
sl = {json.loads(l)["core"]: json.loads(l)["slide"] for l in open(sys.argv[2])}
for r in cen:
    d, s = fit(r); slide = sl[r["core"]]; h, t, _ = batch._headAndTail(slide); hi = max(r["offsets"])
    print(r["core"], slide, "head", h, "tail", t, "fit s", s, "hi", hi, "hi+h-t", hi + h - t, "OK" if s == hi + h - t else "DIFF")
