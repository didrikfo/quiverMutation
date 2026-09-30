"""Overhang fit and the mirror readings of H-021, from census JSONL lines.

Reads the files written by experimentalist_census.py and prints, per core, the
best reflection fit and three readings of "mirror", then the summary tables.

Reflection fit (d = overhang): some s with |s - (lo + hi)| = d, lo/hi the least
and greatest offset with a row, such that every orbit is closed under
o -> s - o wherever s - o is an offset with a row, and at least one orbit holds
some o != s - o together with s - o. Smallest d in 0..6 wins ("no fit" if none).

Mirror readings, for an orbit with held offsets H and mirror offsets M
(M = offsets p whose mirror row lies in the orbit):
  loose   : M nonempty for some orbit.                       (round 001 reading)
  strict  : some orbit has p in M with p not in H, i.e. the mirror of c@p lies
            in the orbit of c@q, and c@p is NOT in that orbit (q != p and the
            placements are in different orbits).            (the STEERING reading)
  strict2 : some orbit has p in M and q in H with q != p.  (weaker: also counts
            a two-offset orbit holding the mirror of one of its own members)

    .venv/bin/python workshop/rounds/002/experimentalist_fit.py logs/t1-census-n13-shard*.jsonl
"""
import json
import sys
from collections import Counter


def fit(record, maxOverhang=6):
    offsets = record["offsets"]
    lo, hi = min(offsets), max(offsets)
    orbits = [set(o["held"]) for o in record["orbits"]]
    for d in range(maxOverhang + 1):
        for s in sorted({lo + hi + d, lo + hi - d}):
            closed, swapped = True, False
            for orbit in orbits:
                for o in orbit:
                    r = s - o
                    if r in offsets:
                        if r not in orbit:
                            closed = False
                        elif r != o:
                            swapped = True
            if closed and swapped:
                return d, s
    return None, None


def readings(record):
    loose = strict = strict2 = False
    for o in record["orbits"]:
        held, mirrors = set(o["held"]), set(o["mirrors"])
        loose |= bool(mirrors)
        strict |= bool(mirrors - held)
        strict2 |= any(p != q for p in mirrors for q in held)
    return loose, strict, strict2


def main(paths):
    records = []
    for path in paths:
        with open(path) as handle:
            records += [json.loads(line) for line in handle if line.strip()]
    records.sort(key=lambda r: r["core"])
    table = Counter()
    for r in records:
        d, s = fit(r)
        loose, strict, strict2 = readings(r)
        small = all(len(o["held"]) <= 2 for o in r["orbits"])
        print("%-8s d=%-4s s=%-4s loose=%d strict=%d strict2=%d small=%d closed=%d orbits=%s" % (
            r["core"], d, s, loose, strict, strict2, small,
            all(o["closed"] for o in r["orbits"]),
            " ".join("{%s}%d" % (",".join(map(str, o["held"])), o["size"]) for o in r["orbits"])))
        table[("fit" if d is not None else "nofit", loose, strict, strict2)] += 1
    print("\ncores: %d; total orbit walks: %d" % (len(records), sum(len(r["orbits"]) for r in records)))
    print("overhang d:", sorted(Counter(fit(r)[0] for r in records).items(), key=lambda x: (x[0] is None, x[0])))
    print("(fit?, loose, strict, strict2) -> cores")
    for k, v in sorted(table.items()):
        print(" ", k, v)


if __name__ == "__main__":
    main(sys.argv[1:])
