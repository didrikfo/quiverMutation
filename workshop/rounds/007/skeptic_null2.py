"""Null for the centre formula s = first + last outside offset (E-060/E-061), n = 13 census + slides.
  .venv/bin/python workshop/rounds/007/skeptic_null2.py [trials]
Null B: shuffle which offsets are 'outside' (same number), keep the fitted s; pred = first+last outside.
Null A: keep the slide, redo the orbit partition at random (same block sizes), refit, compare s to pred.
Reports by class and by informative (>=2 orbits, some pair) vs one-orbit.
"""
import json, os, sys, random, collections
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "workshop/rounds/002"); sys.path.insert(0, "workshop/rounds/007")
from experimentalist_fit import fit
from skeptic_null import nullrec, kind
T = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rng = random.Random(11)
cen = {}; sl = {}
for l in open("workshop/rounds/002/experimentalist_census_n13.jsonl"):
    r = json.loads(l); cen[r["core"]] = r
for l in open("workshop/rounds/003/theorist_slides_n13.jsonl"):
    r = json.loads(l); sl[r["core"]] = r["slide"]
def cls(slide):
    n = len(slide); O = [i for i, x in enumerate(slide) if x == "o"]
    if not O: return "allI"
    if len(O) == n: return "allO"
    a, b = O[0] == 0, O[-1] == n - 1
    return "endboth" if a and b else "end0" if a else "endhi" if b else "int"
res = collections.defaultdict(list)
for c, r in cen.items():
    d, s = fit(r)
    if d is None: continue
    slide = sl[c]; m = len(slide); k = slide.count("o")
    hi = max(r["offsets"]); O = [i for i, x in enumerate(slide) if x == "o"]
    pred = O[0] + O[-1] if O else hi
    # null B: random outside set of same size (offset index = position in slide)
    hitB = 0
    for _ in range(T):
        Oz = sorted(rng.sample(range(m), k)); p = Oz[0] + Oz[-1] if Oz else hi
        hitB += (p == s)
    # null A: random partition, refit, same slide
    hitA = fitA = 0
    for _ in range(T):
        d2, s2 = fit(nullrec(r))
        if d2 is not None:
            fitA += 1; hitA += (s2 == pred)
    res[(cls(slide), kind(r))].append((c, s == pred, hitB / T, fitA / T, hitA / T))
print("class kind n obs_hits  sum P_B(s==pred)  sum P_A(fit&s==pred)  [expected chance hits under each null]")
tot = collections.Counter()
for key, v in sorted(res.items()):
    print(key, len(v), sum(x[1] for x in v), "%.1f" % sum(x[2] for x in v), "%.1f" % sum(x[4] for x in v),
          "fails:", [x[0] for x in v if not x[1]])
print()
for key, v in sorted(res.items()):
    if key[0] == "int": print(key, [(x[0], "%.2f" % x[2], "%.2f" % x[4]) for x in v])
