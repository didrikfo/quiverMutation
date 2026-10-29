"""Round 097 maverick: positive control. Verified-move orbits (lnaMoves.closureUnderMoves, same derived class by an explicit mutation path)
must have constant q_S. Checks the orbit of reps in the split key groups at n = 10, 11, 12.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/097/maverick_orbitcontrol.py N"""
import sys; sys.path.insert(0, "workshop/rounds/097")
import collections
from quivermutation import coxeterTables as ct, lnaMoves
N = int(sys.argv[1]); keys = ct.lnaKeys(N)
sys.argv = ["x", str(N)]
import importlib.util
spec = importlib.util.spec_from_file_location("gs", "workshop/rounds/097/maverick_genusscan.py")
src = open("workshop/rounds/097/maverick_genusscan.py").read().split("split = 0; tot = 0")[0]
g = {"__name__": "gs"}; exec(compile(src, "gs", "exec"), g); inv = g["inv"]
rl0 = next(iter(keys)); print("example key row", rl0, "len", len(rl0))
rep = collections.defaultdict(list)
for rl, k in keys.items(): rep[k].append(rl)
checked = 0; bad = 0; orbsizes = []; evald = 0
for k, v in list(rep.items()):
    for rl in v[:3]:
        orb = lnaMoves.orbitOf(N, list(rl))
        names = orb
        invs = set()
        for nm in names:
            t = tuple(int(c) for c in nm) if isinstance(nm, str) else tuple(nm)
            if t in keys: invs.add(inv(t)); evald += 1
        checked += 1; orbsizes.append(len(names)); bad += len(invs) > 1
print("orbits checked", checked, "orbits with >1 q_S class", bad, "orbit sizes max", max(orbsizes), "members evaluated", evald)
