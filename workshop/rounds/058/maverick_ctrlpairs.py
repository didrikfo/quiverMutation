"""Round 058: pick certified-equivalent control pairs inside class 0 of the n = 10 Phi^18 group (same F-047 move orbit, different rows).
usage (repo root): .venv/bin/python workshop/rounds/058/maverick_ctrlpairs.py"""
import sys, collections
sys.path.insert(0, "workshop/rounds/054"); sys.path.insert(0, "workshop/rounds/029")
import sympy, maverick_pq as mp
from quivermutation import coxeterTables as ct, freeMoves
n = 10
lnas, orbits = freeMoves.derivedOrbits(n, rules=None, free=True, edges=True, doubles=True)
oid = {m: o for o, mem in orbits.items() for m in mem}
r0 = (0,0,0,0,0,0,3,0)
o = orbits[oid[r0]]
print("orbit of 00000030 size", len(o))
for m in sorted(o, key=lambda r:(sum(r), r))[:12]: print("".join(map(str,m)))
