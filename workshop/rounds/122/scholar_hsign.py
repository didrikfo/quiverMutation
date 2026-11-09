"""Round 122 (scholar), referee item 2: the sign of ENS's stabilisation homotopy H (math/0101145v3, proof of Cor. 3.11,
'it is not hard to check', unproved there) in both Leibniz conventions.
S(A) = A * <e1, e2>, d e1 = e2, |e1| = |e2| + 1; A free on x, y, z with d x = y z, d y = d z = 0 (|x| = |y| + |z| + 1).
H(w): take the leftmost ('first') or rightmost ('last') letter of w in {e1, e2}; 0 if it is e1; if it is e2, replace it by e1
with sign (-1)^{|prefix|+1} or (-1)^{|suffix|+1}.  Test i.tau - Id = dH + Hd on every word of length <= 4 (tau kills e1, e2).
Leibniz 'ENS' (left): d(uv) = (du)v + (-1)^|u| u dv.  'engine' (round-106 engine, scholar_ens.py dword): sign from LATER letters.
Usage (repo root): .venv/bin/python -I workshop/rounds/122/scholar_hsign.py"""
from itertools import product

def add(D, w, c):
    D[w] = D.get(w, 0) + c
    if D[w] == 0: del D[w]

def dword(w, d, deg, conv):
    out = {}
    for k, x in enumerate(w):
        s = (-1) ** sum(deg[y] for y in (w[:k] if conv == 'ENS' else w[k + 1:]))
        for u, c in d[x].items(): add(out, w[:k] + u + w[k + 1:], s * c)
    return out

def H(w, deg, pos, sgn):
    idx = [k for k, x in enumerate(w) if x in ('e1', 'e2')]
    if not idx: return {}
    k = idx[0] if pos == 'first' else idx[-1]
    if w[k] == 'e1': return {}
    side = w[:k] if sgn == 'prefix' else w[k + 1:]
    return {w[:k] + ('e1',) + w[k + 1:]: (-1) ** (sum(deg[y] for y in side) + 1)}

def lin(f, E):
    out = {}
    for w, c in E.items():
        for u, e in f(w).items(): add(out, u, c * e)
    return out

fails = {}
for dy, dz, e1 in product([0, 1, 2], [0, 1, 2], [1, 2, 3]):
    deg = {'y': dy, 'z': dz, 'x': dy + dz + 1, 'e1': e1, 'e2': e1 - 1}
    d = {'x': {('y', 'z'): 1}, 'y': {}, 'z': {}, 'e1': {('e2',): 1}, 'e2': {}}
    words = [w for n in range(1, 5) for w in product(['x', 'y', 'z', 'e1', 'e2'], repeat=n)]
    for conv, pos, sgn in product(['ENS', 'engine'], ['first', 'last'], ['prefix', 'suffix']):
        D = lambda w: dword(w, d, deg, conv); h = lambda w: H(w, deg, pos, sgn)
        bad = 0
        for w in words:
            lhs = {} if any(x in ('e1', 'e2') for x in w) else {w: 1}
            add(lhs, w, -1)
            rhs = lin(D, h(w)); [add(rhs, u, c) for u, c in lin(h, D(w)).items()]
            if lhs != rhs: bad += 1
        fails.setdefault((conv, pos, sgn), []).append(bad)
print(f'{len(words)} words x 27 degree choices (|y|,|z| in 0..2, |e1| in 1..3); failures per choice (min..max):')
for k, v in fails.items(): print(f'  Leibniz {k[0]:6} H at {k[1]:5} e-letter, sign (-1)^(|{k[2]}|+1): {min(v)}..{max(v)}' + ('   <- identity holds' if max(v) == 0 else ''))
