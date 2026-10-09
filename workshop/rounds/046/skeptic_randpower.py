"""Round 046 (skeptic): power and soundness check of skeptic_inv.invariants on random unitriangular 7x7 integer matrices (NOT Cartan matrices of algebras).
(a) soundness: C' = P C P^T for random P in GL_7(Z) (product of elementary matrices) keeps every invariant.
(b) power: among random matrices with equal Coxeter charpoly, how often do the other invariants differ?   Usage: skeptic_randpower.py [N]"""
import sys, random
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
sys.path.insert(0, 'workshop/rounds/046'); import skeptic_inv as si
import numpy as np
random.seed(1); n = 7; N = int(ARGV[1]) if len(ARGV) > 1 else 300
def rnd():
    C = np.eye(n, dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.35: C[i, j] = random.choice([1, 1, 1, 2, -1])
    return C
def rndP():
    P = np.eye(n, dtype=int)
    for _ in range(12):
        i, j = random.sample(range(n), 2); E = np.eye(n, dtype=int); E[i, j] = random.choice([1, -1]); P = E @ P
    return P
bad = 0
for _ in range(15):
    C = rnd(); P = rndP(); C2 = P @ C @ P.T
    a, b = si.invariants(C.tolist()), si.invariants(C2.tolist())
    if si.diff(a, b): bad += 1; print('SOUNDNESS FAIL', si.diff(a, b))
print('soundness: 15 congruent pairs, invariant differs in', bad)
groups = {}
for _ in range(N):
    C = rnd(); inv = si.invariants(C.tolist(), mmax=5); groups.setdefault(inv['charpoly'], []).append(inv)
pairs = sep = 0; names = {}
for g in groups.values():
    for a in range(len(g)):
        for b in range(a + 1, len(g)):
            pairs += 1; d = si.diff(g[a], g[b])
            if d: sep += 1
            for k in d: names[k] = names.get(k, 0) + 1
print('random matrices', N, 'charpoly groups', len(groups), 'pairs with equal charpoly', pairs, 'separated by another invariant', sep)
print({k: v for k, v in names.items()})
