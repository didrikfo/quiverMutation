"""Round 086 (idea): within each key group at n = N (N=12) with indefinite periodic members: invariants (eigenvalue-signature profile, SNF(C+C^T)-type
Cartan data), then full-lattice equivariant-isometry (r082 maverick_fulliso.isometryH) between one rep per (profile, invariant) cell and a
second member of the same cell/different cell.  A NO with node cap not hit = Cartan-incongruent (proof for that pair) with equal key and profile.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/086/maverick_pairs.py N [--plan]"""
import sys, time, collections, pickle, random
sys.path.insert(0, "workshop/rounds/074"); sys.path.insert(0, "workshop/rounds/082")
import numpy as np
from maverick_isometry import data
from maverick_eigsig import sigProfile
from maverick_fulliso import majorant, isometryH
import os, tempfile
N = int(sys.argv[1]); t0 = time.time()
d = pickle.load(open(os.environ.get("MAVERICK_PKL", os.path.join(tempfile.gettempdir(), "maverick_scan%d.pkl" % N)), "rb")); multi, per = d["multi"], d["per"]
cells = collections.defaultdict(list)   # (key index, profile) -> members, only indefinite periodic
for ki, (k, v) in enumerate(multi.items()):
    for rl in v:
        neg, p = per[rl]
        if neg > 0 and p > 0:
            G, M, _ = data(rl, N); cells[(ki, sigProfile(G, M))].append(rl)
print("cells (key,profile) with >=2 indefinite periodic members:", sum(len(v) > 1 for v in cells.values()), "of", len(cells), "sizes", sorted(len(v) for v in cells.values())[-10:], flush=True)
if "--plan" in sys.argv: sys.exit()
import signal
class TO(Exception): pass
def h(sig, fr): raise TO()
signal.signal(signal.SIGALRM, h)
K = 10**9 if 'all' in sys.argv else (int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 6); LIM = 90; tally = collections.Counter()
for (ki, prof), mem in sorted(cells.items(), key=lambda kv: len(kv[1])):
    neg = {float(a) for a, pp, q in prof if q}; assert all(pp == 0 for a, pp, q in prof if q)
    random.seed(7); others = random.sample(mem[1:], min(K, len(mem) - 1)); a = mem[0]
    Ga, Ma, _ = data(a, N); Ha = majorant(Ga, Ma, neg)
    for b in others:
        t1 = time.time(); Gb, Mb, _ = data(b, N)
        try:
            signal.alarm(LIM); Hb = majorant(Gb, Mb, neg); res, P, nodes, vs = isometryH(Ga, Ma, Gb, Mb, Ha, Hb); signal.alarm(0)
        except TO: res, nodes, vs = "TIMEOUT%ds" % LIM, 0, []
        tally[(ki, res)] += 1
        if res != "YES" or K < 100: print("group", ki, "".join(map(str, a)), "->", "".join(map(str, b)), res, "gens", len(vs), "nodes", nodes, f"{time.time()-t1:.0f}s", flush=True)
print("tally (group, result):", dict(tally), f"{time.time()-t0:.0f}s")
