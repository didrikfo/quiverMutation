"""Round 071 (idea): number of orbit-plus-mirror classes (freeMoves.derivedOrbits(free=True) plus mirrorRow union) inside each
n = 11 periodic-indefinite group (E-190). Unions are certified derived equivalences, so #classes in a group is an upper bound
on its number of derived classes (as far as the moves reach). Also lists where the 16 unplaced (E-162) LNAs lie.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/071/maverick_orbits.py [N]"""
import sys, time, collections
import numpy as np
from quivermutation import coxeterTables as ct, nakayama as nk, freeMoves as fm
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
t0 = time.time()
lnas, orbits = fm.derivedOrbits(N, free=True)
print(f"derivedOrbits {len(orbits)} orbits of {len(lnas)}  {time.time()-t0:.0f}s", flush=True)
parent = {}
for k, mem in orbits.items():
    for m in mem: parent[tuple(m)] = tuple(mem[0])
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for m in list(parent):
    a, b = find(m), find(fm.mirrorRow(N, m))
    if a != b: parent[a] = b
def periodic(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return True
        M = M @ Phi
    return False
keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl in lnas:
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    if not periodic(np.rint(-Ci.T @ C).astype(np.int64)): continue
    groups[keys[rl]].append(tuple(rl))
size = collections.Counter(find(m) for m in parent)
for key, mem in sorted(groups.items(), key=lambda kv: len(kv[1])):
    cls = collections.Counter(find(m) for m in mem)
    whole = [size[c] for c in cls]
    print(f"group {len(mem)}: orbit+mirror classes {len(cls)}; members per class {sorted(cls.values(), reverse=True)[:12]}; whole-class sizes {sorted(whole, reverse=True)[:12]}")
    for r, c in sorted(cls.items(), key=lambda kv: -kv[1])[:8]: print("   rep", ct.className(r), "members in group", c, "class size", size[r])
print(f"total {time.time()-t0:.0f}s")
