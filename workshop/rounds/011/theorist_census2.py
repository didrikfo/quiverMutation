"""Census of the --max-word 4 catalogue at length n using only the orbits of the SINGLE-RELATION rows k@o.
Walk every single-relation row's reduced orbit (k = 3..n-2, all offsets).  A catalogue word is 'predicted parity-class' when
every placement's reduced row lies in a single-relation orbit and the orbit ids at consecutive offsets alternate a, b, a, b
(two orbits, offset parity).  Words with a placement in no single-relation orbit are 'unknown' (need their own walk).
The two orbits must carry one Coxeter key and each hold its own mirror (406 fails the last: its two orbits are mirror images).
Compare with the E-076 lists A (even n) / B (odd n).
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_census2.py N [limit]
"""
import sys, time
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm, coxeterTables as ct
A = "35 455 3334 3336 5003 5055 5504 5505 5506".split()
B = "36 405 466 3335 5004 5006 5046 5056 5066 5605".split()
n = int(sys.argv[1]); limit = int(sys.argv[2]) if len(sys.argv) > 2 else 400000
t = time.time(); oid = {}; sizes = []; capped = 0; okey = {}; omir = {}
for k in range(3, n - 1):
    for o in range(n):
        row = batch._rowFor(n, str(k), o)
        if row is None: continue
        s = fm._startOf(tuple(row), fm.REDUCED)
        if s in oid: continue
        rep = fm.orbitReport(n, tuple(row), free=fm.REDUCED, limit=limit)
        if not rep.closed: capped += 1; continue
        i = len(sizes); sizes.append(len(rep.rows)); okey[i] = ct.lnaCoxeterKey(n, tuple(row))
        omirrow = fm._startOf(fm.mirrorRow(n, tuple(row)), fm.REDUCED); omir[i] = omirrow in rep.rows
        for r in rep.rows: oid[r] = i
print('n = %d: %d single-relation orbits walked, %d capped, %.0fs' % (n, len(sizes), capped, time.time() - t), flush=True)
pred, unknown, other = [], 0, []
for w in batch._singleCores(4, 6, False):
    ids = []
    for o in range(n):
        row = batch._rowFor(n, w, o)
        if row is None: continue
        ids.append(oid.get(fm._startOf(tuple(row), fm.REDUCED)))
    if len(ids) < 2: continue
    if None in ids: unknown += 1; continue
    if len(set(ids)) == 2 and len({okey[i] for i in ids}) == 1 and all(omir[i] for i in ids) and all(ids[i] != ids[i + 1] for i in range(len(ids) - 1)): pred.append(w)
    elif len(set(ids)) > 1 and len(set(ids)) <= 3 : other.append((w, ids))
ref = A if n % 2 == 0 else B
print('predicted (alternating between two single-relation orbits):', ' '.join(pred))
print('E-076 list:', ' '.join(ref))
print('in prediction not in list:', sorted(set(pred) - set(ref)), '; in list not in prediction:', sorted(set(ref) - set(pred)))
print('words with a placement in no single-relation orbit: %d' % unknown)
print('words in 2..3 single-relation orbits, not alternating:', [(w, ''.join(map(str, [chr(48 + x % 10) for x in ids]))) for w, ids in other][:12], len(other))
