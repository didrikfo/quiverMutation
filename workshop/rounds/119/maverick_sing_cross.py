"""Round 119 maverick: more cross-class pairs of the singular n = 12 group of 47 (SNF(S) = (1^10, 20) for 14 members vs (1^10, 12) for 33): E-269-style commutant search must find no P.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_sing_cross.py [number of pairs = 3] [B = 100]"""
import sys, time, collections; sys.path.insert(0, "workshop/rounds/119"); sys.path.insert(0, "workshop/rounds/097")
import maverick_sing as W
from maverick_cartanZ import cartan
from quivermutation import coxeterTables as ct
keys = ct.lnaKeys(12); kk = keys[(2, 2, 4, 0, 0, 3, 3, 0, 3, 0)]; mem = [tuple(r) for r, k in keys.items() if k == kk]
fp = collections.defaultdict(list)
for rl in mem: fp[W.qbar_fp(cartan(12, list(rl)))].append(rl)
A, B = sorted(fp.values(), key=len)
npairs = int(sys.argv[1]) if len(sys.argv) > 1 else 3; BB = int(sys.argv[2]) if len(sys.argv) > 2 else 100
for i in range(npairs):
    a, b = A[(3 * i) % len(A)], B[(5 * i + 1) % len(B)]; t = time.time()
    r = W.search(cartan(12, list(a)), cartan(12, list(b)), [BB])
    print("cross", a, b, [{k: v for k, v in x.items() if k != "P"} for x in r], "%.0fs" % (time.time() - t), flush=True)
