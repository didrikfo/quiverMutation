"""For a word at n: does each parity orbit hold the mirror of its own placements, or of the other's?
  .venv/bin/python workshop/rounds/011/theorist_mirror406.py N word"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm
n = int(sys.argv[1]); w = sys.argv[2]
orbs = []
for o in (0, 1):
    rep = fm.orbitReport(n, tuple(batch._rowFor(n, w, o)), free=fm.REDUCED, limit=1500000); orbs.append(rep.rows)
for o in range(n):
    row = batch._rowFor(n, w, o)
    if row is None: continue
    m = fm._startOf(fm.mirrorRow(n, tuple(row)), fm.REDUCED)
    print(w, o, 'in orbit of offset', o % 2, ': own', fm._startOf(tuple(row), fm.REDUCED) in orbs[o % 2], '; mirror in orbit0', m in orbs[0], 'orbit1', m in orbs[1])
