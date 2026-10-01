"""All single-relation rows k@o at length n grouped by Coxeter key: orbits in each key, polynomial coefficients.
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_keyclass.py N word [word..]
Prints for the key of each given word@offset: the key, the single-relation rows in that key class,
their reduced-walk orbit sizes."""
import sys
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm, coxeterTables as ct
n = int(sys.argv[1])
rows = {}
for k in range(2, n - 1):
    for o in range(n):
        r = batch._rowFor(n, str(k), o)
        if r is not None: rows[(k, o)] = tuple(r)
bykey = {}
for (k, o), r in rows.items(): bykey.setdefault(ct.lnaCoxeterKey(n, r), []).append((k, o))
for w in sys.argv[2:]:
    for o in range(n):
        r = batch._rowFor(n, w, o)
        if r is None: continue
        key = ct.lnaCoxeterKey(n, tuple(r))
        print(w, '@', o, 'key', key)
        sizes = {}
        for (k, oo) in bykey.get(key, []):
            s = fm._startOf(rows[(k, oo)], fm.REDUCED)
            rep = fm.orbitReport(n, rows[(k, oo)], free=fm.REDUCED, limit=400000)
            sizes[(k, oo)] = len(rep.rows)
        print('   single-relation rows with this key and their orbit sizes:', sizes)
        break
