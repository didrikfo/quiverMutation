"""Round 102 (toolsmith, response to referee item 2): give isoTo (Sect. 6 target test) a NO on KEY-EQUAL non-isomorphic targets.
usage: .venv/bin/python workshop/rounds/102/toolsmith_isoNO.py"""
import sys, warnings; warnings.simplefilter('ignore')
sys.path.insert(0, 'workshop/rounds/102'); sys.path.insert(0, '.'); sys.path.insert(0, 'tests')
import toolsmith_litex as L
from helpers import path_algebra, quiet
for n, first in ((3, [1, 3]), (4, [1, 3, 5])):
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()): c = L.sect6(n, first)
    base = L.key(L.a2n3(n)); R = L.rect(n)
    print('n =', n, 'key(rect) == key(start):', L.key(R) == base)
    cands = {'pushed-square quiver (after steps %s)' % first: c, 'start A_{2n,3}': L.a2n3(n)}
    # rectangle arrows with the first commutativity relation replaced by a zero relation on one path / dropped
    ar = sorted(L.arrow_set(R)); rl = [list(r) for r in sorted(L.rel_set(R))]
    cands['rect, relations dropped'] = path_algebra(ar, [])
    cands['rect, first square zero on one path'] = path_algebra(ar, [[rl[0][0]]] + rl[1:])
    cands['rect, all squares both paths zero'] = path_algebra(ar, [[p] for r in rl for p in r])
    for name, t in cands.items():
        k = L.key(t)
        print('  %-45s key equal: %-5s iso to rect: %s' % (name, k == base, L.isoTo(t, R) is not None))
