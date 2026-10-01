"""Letter statistics of the reduced orbits P (3@2 / 3@3) and Q (5@0 / 6@0).  .venv/bin/python workshop/rounds/013/theorist_orbitstats.py N"""
import sys
sys.path.insert(0, '.')
import batch
from collections import Counter
from quivermutation import freeMoves as fm
n = int(sys.argv[1])
for name, (k, o) in (('P', ('3', 2 if n % 2 == 0 else 3)), ('Q', ('5' if n % 2 == 0 else '6', 0))):
    rep = fm.orbitReport(n, tuple(batch._rowFor(n, k, o)), free=fm.REDUCED, limit=300000)
    rows = rep.rows
    nrel = Counter(sum(1 for x in r if x) for r in rows)
    mx = Counter(max(r) for r in rows)
    mn = Counter(min([x for x in r if x] or [0]) for r in rows)
    ssum = Counter(sum(r) for r in rows)
    print(n, name, len(rows), 'nrel', sorted(nrel.items()), 'max', sorted(mx.items()), 'min', sorted(mn.items()))
    print('   sum of letters', sorted(ssum.items()))
    print('   rows with a single relation:', sorted(''.join(map(str, r)) for r in rows if sum(1 for x in r if x) == 1))
