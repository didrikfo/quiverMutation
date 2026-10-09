"""T6 (round 023, maverick): LNAs with >= 2 relations of >= 3 arrows at n = 8, 9: D1 and peeling formula against the recorded
minimum cord depth (workshop/rounds/022/theorist_blocked_depths.txt: all LNAs of depth >= 2, n = 8..10).
Depth 1 for every LNA not listed there (D1 is E-103's, 0 mismatches over all LNAs, so unlisted = depth 1 iff some big relation unblocked).
usage: maverick_twobig.py N   (repository root)"""
import sys, itertools
sys.path.insert(0, 'workshop/rounds/022')
N = int(sys.argv[1])
data = {}
for line in open('workshop/rounds/022/theorist_blocked_depths.txt'):
    n, d, _, dep = line.split()
    data[(int(n), d)] = int(dep)
def formula(d):
    r = [0] + list(map(int, d)) + [0] * 5; best = None
    for s in range(1, len(d) + 1):
        m = r[s]
        if m < 3: continue
        a = 0
        while s - 1 - a >= 1 and r[s - 1 - a] == 2: a += 1
        b = 0
        while r[s + m - 1 + b] == 2: b += 1
        v = 1 + min(a, b); best = v if best is None else min(best, v)
    return best
def nbig(d): return sum(1 for c in d if int(c) >= 3)
from quivermutation import coxeterTables as ct
lst = ["".join(map(str, r)) for r in sorted(ct.lnaStatus(N))]
two = [d for d in lst if nbig(d) >= 2]
print('n', N, 'LNAs', len(lst), 'with >=2 big', len(two), 'of which listed depth>=2', sum((N, d) in data for d in two))
for d in two:
    dep = data.get((N, d), 1); f = formula(d)
    print(d, 'depth', dep, 'formula', f, 'OK' if f == dep else 'FAIL')
