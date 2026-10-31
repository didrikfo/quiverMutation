"""Round 103 literature check: compare the move orbit of B_n (E-247) with the published table of arXiv:2305.06642
(section 'Complete classification for n <= 8', transcribed in tests/paper_classification.py), and check the
published A_n ~ B_n (Cor. 4.6 line 3, r = 0) against the repo's move orbit.
Row convention (repo): n-2 entries, entry i = number of arrows of the relation starting at vertex i+1 (0 = none).
usage: .venv/bin/python workshop/rounds/103/experimentalist_fosse_check.py"""
import sys
sys.path.insert(0, 'tests')
from paper_classification import PAPER_CLASSES
from quivermutation import freeMoves as fm
A = lambda n: tuple([n - 1] + [0] * (n - 3))
B = lambda n: tuple([0] * (n - 4) + [3, 0])
def rels(row):
    return [(i + 1, l) for i, l in enumerate(row) if l]
def sep3(row):
    """relations of length >= 3 only, almost separate (consecutive ones share <= 1 vertex); None if not."""
    r = [(s, l) for s, l in rels(row) if l >= 3]
    for (s1, l1), (s2, l2) in zip(r, r[1:]):
        if s2 < s1 + l1:   # shares >= 2 vertices (end vertex of first is s1+l1)
            return None
    return tuple(r)
for n in range(6, 9):
    cls = {}
    for name, lst in PAPER_CLASSES[n].items():
        for (st, ln) in lst:
            cls[tuple(zip(st, ln))] = name
    orb = fm.orbitReport(n, B(n), free=True, limit=200000).rows
    seen = {}
    for row in orb:
        s = sep3(row)
        if s is None: continue
        seen.setdefault(s, 0); seen[s] += 1
    names = {cls.get(s, 'NOT IN TABLE') for s in seen}
    print('n=%d |orbit B_n|=%d  distinct almost-separate (>=3) shapes in orbit: %d  table classes hit: %s' % (n, len(orb), len(seen), sorted(names)))
    print('   shapes:', sorted(seen))
    dn = [k for k, v in cls.items() if v.startswith('D_') ]
    print('   table D_%d row:' % n, sorted(dn), ' A_n shape', ((1, n - 1),), 'in table D class:', cls.get(((1, n - 1),)), ' A_n in orbit(B_n):', A(n) in orb)
