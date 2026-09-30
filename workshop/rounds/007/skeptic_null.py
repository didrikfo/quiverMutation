"""Round 007 null test for |R|<=4 reflection fits (E-060/E-061/E-062 data).
  .venv/bin/python workshop/rounds/007/skeptic_null.py [trials]
Null: keep each core's offsets and orbit block sizes, randomly reassign offsets to blocks
(uniform over labelled permutations), run the project's own fit() (d<=6).
Only the plain uniform null is implemented (a neighbour-aware null is not).
"""
import json, os, sys, random, collections, glob
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "workshop/rounds/002")
from experimentalist_fit import fit
T = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rng = random.Random(7)

def nullrec(r):
    offs = list(r["offsets"]); rng.shuffle(offs)
    out, i = [], 0
    for o in r["orbits"]:
        k = len(o["held"]); out.append({"held": offs[i:i+k]}); i += k
    return {"offsets": r["offsets"], "orbits": out}

def kind(r):
    sizes = sorted(len(o["held"]) for o in r["orbits"])
    if len(sizes) == 1: return "one-orbit"
    if max(sizes) == 1: return "all-singleton"
    return "informative"

def run(path, label):
    rows = []
    for l in open(path):
        r = json.loads(l)
        # held must partition offsets
        if sorted(x for o in r["orbits"] for x in o["held"]) != sorted(r["offsets"]): 
            cov = "nonpartition"
        else: cov = "ok"
        d, s = fit(r)
        if d is None: rows.append((r["core"], kind(r), None, None, None, None, cov)); continue
        pf = pd = 0
        for _ in range(T):
            d2, s2 = fit(nullrec(r))
            if d2 is not None:
                pf += 1
                if d2 <= d: pd += 1
        rows.append((r["core"], kind(r), d, s, pf / T, pd / T, cov))
    return rows

if __name__ == "__main__":
    files = [("n13", "workshop/rounds/002/experimentalist_census_n13.jsonl")]
    files += [("n15", "workshop/rounds/004/experimentalist_census_12cores_n15.jsonl"),
              ("n16", "workshop/rounds/004/experimentalist_census_12cores_n16.jsonl")]
    for lab, p in files:
        rows = run(p, lab)
        print("==", lab, "cores", len(rows), "T", T)
        print("core kind d s P(null fit) P(null d<=d_obs) partition")
        for x in rows: print(*x)
        fits = [x for x in rows if x[2] is not None]
        print("fits:", len(fits), collections.Counter(x[1] for x in fits))
        for k in ("one-orbit", "informative"):
            v = [x for x in fits if x[1] == k]
            if not v: continue
            for a in (0.05, 0.2):
                print(k, "n=%d" % len(v), "p(d<=dobs) < %.2f: %d" % (a, sum(x[5] < a for x in v)),
                      "mean P(null fit) %.3f" % (sum(x[4] for x in v) / len(v)),
                      "expected chance fits %.1f" % sum(x[4] for x in v))
