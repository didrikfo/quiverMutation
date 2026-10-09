"""T6 (round 023, maverick): distribution of the mirror-chain depth over two-big-relation LNAs, n = 8..12, and the ones with predicted depth >= 3.
usage: maverick_predict2.py   (repository root)"""
import sys, collections
sys.path.insert(0, 'workshop/rounds/023')
from maverick_chain import formula
from quivermutation import coxeterTables as ct
for n in (8, 9, 10, 11, 12):
    lst = ["".join(map(str, r)) for r in sorted(ct.lnaStatus(n))] if n <= 11 else []
    c = collections.Counter(); ex = []
    for d in lst:
        if sum(1 for ch in d if int(ch) >= 3) < 2: continue
        v = formula(d, lambda l: True); c[v] += 1
        if v >= 3: ex.append(d)
    print(n, dict(c), ex[:10])
