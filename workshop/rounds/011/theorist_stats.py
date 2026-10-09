"""Integer statistics of rows over a reduced-walk orbit: which are constant (mod 2) on an orbit?
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_stats.py N word [word..]
"""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm

def stats(row):
    n = len(row) + 2
    rel = [(i + 1, i + 1 + r) for i, r in enumerate(row) if r]
    return {
        'count': len(rel),
        'sumlen': sum(r for r in row),
        'sumstart': sum(a for a, b in rel),
        'sumend': sum(b for a, b in rel),
        'sum(start+end)': sum(a + b for a, b in rel),
        'sum_start*len': sum(a * (b - a) for a, b in rel),
        'nOdd': sum(1 for r in row if r % 2),
        'sumOddPos': sum(i + 1 for i, r in enumerate(row) if r % 2),
        'sumEvenPos': sum(i + 1 for i, r in enumerate(row) if r and r % 2 == 0),
    }

for n in [int(sys.argv[1])]:
    for w in sys.argv[2:]:
        for o in range(n):
            row = batch._rowFor(n, w, o)
            if row is None: continue
            rep = fm.orbitReport(n, tuple(row), free=fm.REDUCED, limit=400000)
            vals = {}
            for r in rep.rows:
                for k, v in stats(r).items(): vals.setdefault(k, set()).add(v)
            out = []
            for k, s in vals.items():
                for mod in (2, 4):
                    c = {v % mod for v in s}
                    if len(c) == 1: out.append('%s mod %d = %d' % (k, mod, c.pop()))
            print(w, o, len(rep.rows), '; '.join(out), flush=True)
