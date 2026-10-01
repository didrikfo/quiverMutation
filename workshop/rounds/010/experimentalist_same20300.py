"""n = 16: is the size-20300 orbit of 348 {2},{3} and 349 {1},{3} the 4056 orbit (E-064)?
Walks each offset's reduced orbit (limit 1500000), prints size, closed, and pairwise shared rows,
plus whether a start row (or its mirror) of one lies in the other.
  timeout 10m .venv/bin/python workshop/rounds/010/experimentalist_same20300.py
"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm

n = 16
want = [('4056', 1), ('4056', 2), ('348', 2), ('348', 3), ('349', 1), ('349', 3)]
walks = {}
for w, o in want:
    row = tuple(batch._rowFor(n, w, o))
    rep = fm.orbitReport(n, row, free=fm.REDUCED, limit=1500000)
    walks[(w, o)] = (set(rep.rows), fm._startOf(row, fm.REDUCED), fm._startOf(fm.mirrorRow(n, row), fm.REDUCED))
    print(w, o, len(rep.rows), rep.closed, flush=True)
keys = list(walks)
for i, a in enumerate(keys):
    for b in keys[i + 1:]:
        A, sa, ma = walks[a]; B, sb, mb = walks[b]
        print("%s vs %s: shared %d; a-start in b %s; a-mirror in b %s" % (a, b, len(A & B), sa in B, ma in B), flush=True)
