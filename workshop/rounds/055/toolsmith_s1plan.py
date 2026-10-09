"""S-1 n = 15 sizing (round 055). usage: python workshop/rounds/055/toolsmith_s1plan.py N [SAMPLE]
Times enumeration of allRelationLengths(N) and lnaCoxeterKey per row; extrapolates row count and scan time."""
import sys, time
from quivermutation import coxeterTables as ct, nakayama
n = int(sys.argv[1]); S = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
t = time.time(); rows = []
for i, r in enumerate(nakayama.allRelationLengths(n)):
    if i >= S: break
    rows.append(tuple(r))
te = time.time() - t
t = time.time()
for r in rows[:5000]: ct.lnaCoxeterKey(n, r)
tk = (time.time() - t) / min(5000, len(rows))
print('n', n, 'enum %.2f us/row' % (te / len(rows) * 1e6), 'key %.2f ms/row' % (tk * 1e3))
