"""Round 058: list the n = 10 Phi^18 group (E-172): classes, sizes, SNF profile sets, short representatives, and pairs with disjoint profile sets.
usage (repo root): .venv/bin/python workshop/rounds/058/maverick_group.py"""
import sys, collections
sys.path.insert(0, "workshop/rounds/054"); sys.path.insert(0, "workshop/rounds/029")
import sympy, maverick_pq as mp, toolsmith_snfresolve as sr
from quivermutation import coxeterTables as ct
n = 10; keys, comps, *_ = mp.classesOf(n)
by = collections.defaultdict(list)
for c, mem in comps.items(): by[keys[mem[0]]].append(mem)
for k, cl in by.items():
    if len(cl) < 2: continue
    C = sympy.Matrix(ct.lnaCartanMatrix(n, cl[0][0])); Phi = -(C.inv().T) * C
    if (Phi**18) != sympy.eye(n): continue
    sets = [{sr.profile(n, m) for m in mem} for mem in cl]
    for i, mem in enumerate(cl):
        reps = sorted(mem, key=lambda r: (sum(r), r))[:4]
        print("class", i, "size", len(mem), "profiles", len(sets[i]), "reps", ["".join(map(str, r)) for r in reps])
    for i in range(len(cl)):
        for j in range(i + 1, len(cl)):
            print("pair", i, j, "profile sets disjoint:", not (sets[i] & sets[j]), "overlap", len(sets[i] & sets[j]))
