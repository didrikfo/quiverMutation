"""Which GF(2)-affine functionals of a row are constant on a reduced-walk orbit? (T1/T3, E-076)
  timeout 10m .venv/bin/python workshop/rounds/011/theorist_invariant.py N word [word..]
Features per row: r_i mod 2 and [r_i != 0] for i = 1..n-2, and i*r_i etc. are linear combos of these.
Prints, for each placement, orbit size, and the dimension of the invariant space, and whether the
invariant space separates offsets (i.e. the offsets of the word fall into >1 class).
"""
import sys, itertools
sys.path.insert(0, '.')
import batch
from quivermutation import freeMoves as fm

def feat(row):
    return [r & 1 for r in row] + [1 if r else 0 for r in row]

def nullInvariants(rows):
    """return basis of functionals f (vectors over GF(2)) with f.(x-x0)=0 for all x in rows."""
    rows = list(rows); x0 = feat(rows[0]); m = len(x0)
    diffs = []
    for r in rows[1:]:
        d = [a ^ b for a, b in zip(feat(r), x0)]
        diffs.append(int(''.join(map(str, d)), 2))
    # row-reduce diffs, then nullspace
    basis = {}
    for d in diffs:
        while d:
            h = d.bit_length() - 1
            if h in basis: d ^= basis[h]
            else: basis[h] = d; break
    piv = set(basis)
    # nullspace of the matrix with rows = basis vectors; reduce fully
    keys = sorted(basis, reverse=True)
    for k in keys:
        for k2 in keys:
            if k2 != k and (basis[k2] >> k) & 1: basis[k2] ^= basis[k]
    free = [c for c in range(m) if c not in piv]
    null = []
    for f in free:
        v = 1 << f
        for k in keys:
            if (basis[k] >> f) & 1: v |= 1 << k
        null.append(v)
    return null, m

def pf(v, m):
    names = ['p%d' % (i + 1) for i in range(m // 2)] + ['z%d' % (i + 1) for i in range(m // 2)]
    return '+'.join(names[m - 1 - b] for b in range(m) if (v >> b) & 1)

if __name__ == '__main__':
    n = int(sys.argv[1])
    for w in sys.argv[2:]:
        for o in range(n):
            row = batch._rowFor(n, w, o)
            if row is None: continue
            rep = fm.orbitReport(n, tuple(row), free=fm.REDUCED, limit=400000)
            null, m = nullInvariants(rep.rows)
            print(w, o, 'size', len(rep.rows), 'closed', rep.closed, 'dim', len(null), flush=True)
            if len(null) <= 6:
                for v in null: print('    ', pf(v, m))
