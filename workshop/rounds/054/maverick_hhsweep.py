"""Hochschild cohomology dims of every LNA, n = 3..N (round 054), and a non-tree positive control for the code.
usage: python workshop/rounds/054/maverick_hhsweep.py N"""
import sys, collections
sys.path.insert(0, "workshop/rounds/054")
import maverick_pq as mp
from quivermutation import nakayama
# positive control: incidence algebras of posets have HH^m = H^m(order complex).  crown a,b < c,d (circle) -> (1,1);
# the chain a<b<c plus a<d<c (diamond-like, contractible) -> (1).
def inc(rel, n): return {(a, a) for a in range(1, n + 1)} | set(rel)
print("crown (circle)   ", mp.hh(4, None, inc([(1, 3), (1, 4), (2, 3), (2, 4)], 4)), "expect (1, 1)")
print("commutative sq   ", mp.hh(4, None, inc([(1, 2), (1, 3), (1, 4), (2, 4), (3, 4)], 4)), "expect (1,)")
# 2 minimal, 3 maximal, complete bipartite: order complex is K_{2,3}, b1 = 2
print("K_{2,3}          ", mp.hh(5, None, inc([(a, b) for a in (1, 2) for b in (3, 4, 5)], 5)), "expect (1, 2)")
for n in range(3, int(sys.argv[1]) + 1):
    c = collections.Counter(mp.hh(n, rl) for rl in nakayama.allRelationLengths(n))
    print("n =", n, "LNAs", sum(c.values()), "HH profiles", {tuple(int(x) for x in k): v for k, v in c.items()})
