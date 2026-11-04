"""Round 115 (theorist): planted-NO controls for the HH invariant.  Pairs with the SAME Cartan matrix (hence same Coxeter key and same Cartan-level
invariants) that Happel's theorem makes derived-inequivalent exactly when HH^m differ.  Also positive controls (known derived-equivalent: HH must agree)."""
import sys, numpy as np
sys.path.insert(0, '.')
exec(open('workshop/rounds/115/theorist_hh.py').read().split("if __name__")[0])
def cartan(V, edges, rels):
    A = Alg(V, edges, rels); ix = {v: i for i, v in enumerate(A.V)}; C = np.zeros((len(V), len(V)), dtype=int)
    for p, s, t in A.basis: C[ix[s], ix[t]] += 1
    return C, A.n
sq = [(1, 2, 0), (1, 3, 0), (2, 4, 0), (3, 4, 0)]; ab = ((1, 2, 0), (2, 4, 0)); cd = ((1, 3, 0), (3, 4, 0))
def tail(k): return [(4 + i, 5 + i, 0) for i in range(k)]
cases = []
for k in (0, 3):   # n = 4 and n = 7
    V = list(range(1, 5 + k)); E = sq + tail(k)
    cases.append(('comm square +tail%d' % k, V, E, [{ab: 1, cd: -1}]))
    cases.append(('square ab=0 +tail%d' % k, V, E, [{ab: 1}]))
    cases.append(('comm square (+) +tail%d, coeff +' % k, V, E, [{ab: 1, cd: 1}]))
res = {}
for name, V, E, R in cases:
    C, d = cartan(V, E, R); h = hh(V, E, R, top=2)
    res[name] = (C, h); print(name, 'dim', d, 'cartan-sum', C.sum(), 'HH0-2', h[1])
for k in (0, 3):
    a = res['comm square +tail%d' % k]; b = res['square ab=0 +tail%d' % k]
    print('tail', k, 'Cartan equal:', (a[0] == b[0]).all(), ' HH equal:', a[1][1] == b[1][1], '-> HH says', 'YES' if a[1][1] == b[1][1] else 'NO')
# positive controls: comm square (-) vs (+) are isomorphic (sign change on one arrow) -> HH equal; D4 hereditary vs comm square derived-equivalent
print('D4 (hereditary)', hh([1, 2, 3, 4], [(1, 2, 0), (1, 3, 0), (1, 4, 0)], [], top=2)[1])
print('A4 hereditary', hh([1, 2, 3, 4], [(1, 2, 0), (2, 3, 0), (3, 4, 0)], [], top=2)[1])
print('Kronecker x A2-chain', hh([1, 2, 3], [(1, 2, 0), (1, 2, 1), (2, 3, 0)], [], top=2)[1])
