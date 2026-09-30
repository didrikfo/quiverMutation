"""Round 004 extension of rounds/003/theorist_shortfall.py: the interior / end-touch split as a column.

  .venv/bin/python workshop/rounds/004/theorist_shortfall.py [--all]
Reads the round-002 n=13 census (orbits) and round-003 slides. For each core with a reflection fit:
  cls   : allI (no outside offset) | allO (every offset outside) | int (outside block touches neither end) |
          end0 / endhi / endboth (block touches offset 0 / hi / both -- but not allO)
  pred  : first outside + last outside offset (hi if no outside offset)
  s     : fitted centre (experimentalist_fit.fit); ok = (s == pred)
  cons  : slide-consistent centres s' (slide[o]==slide[s'-o] wherever both in R, some nontrivial pair), s' in 0..2hi
  in    : is the fitted s in cons
"""
import json, os, sys, collections
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "workshop/rounds/002")
from experimentalist_fit import fit
cen = {}
for l in open("workshop/rounds/002/experimentalist_census_n13.jsonl"):
    r = json.loads(l); cen[r["core"]] = r
sl = {}
for l in open("workshop/rounds/003/theorist_slides_n13.jsonl"):
    r = json.loads(l); sl[r["core"]] = r["slide"]

def cls(slide):
    n = len(slide); O = [i for i, x in enumerate(slide) if x == "o"]
    if not O: return "allI"
    if len(O) == n: return "allO"
    a, b = O[0] == 0, O[-1] == n - 1
    return "endboth" if a and b else "end0" if a else "endhi" if b else "int"

def consistent(slide):
    n = len(slide); out = []
    for s in range(0, 2 * n - 1):
        pairs = [(o, s - o) for o in range(n) if 0 <= s - o < n]
        if all(slide[o] == slide[r] for o, r in pairs) and any(o != r for o, r in pairs): out.append(s)
    return out

rows = []
for c, r in cen.items():
    d, s = fit(r)
    if d is None: continue
    slide = sl[c]; O = [i for i, x in enumerate(slide) if x == "o"]
    hi = max(r["offsets"]); pred = O[0] + O[-1] if O else hi
    cons = consistent(slide)
    rows.append((c, slide, cls(slide), s, hi, hi - s, pred, s == pred, s in cons, cons))
if "--all" in sys.argv or True:
    print("core slide cls s hi d pred ok in_cons cons")
    for x in rows: print(*x[:9], ",".join(map(str, x[9])))
tab = collections.defaultdict(list)
for x in rows: tab[x[2]].append(x)
print()
for k, v in sorted(tab.items()):
    print(k, "n=%d" % len(v), "pred ok %d" % sum(x[7] for x in v), "s in cons %d" % sum(x[8] for x in v),
          "pred fails:", [x[0] for x in v if not x[7]])
