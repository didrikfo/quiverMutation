"""Round 017 experimentalist: run workshop/rounds/014/scholar_walk.py (E-084 guarded walk) with
arrowPaths.reduceAgainstPivots replaced by the FULL reduction of workshop/rounds/015/theorist_fix.py
(library untouched). Same arguments as scholar_walk.py; `--unpatched` runs it unchanged (control)."""
import sys
sys.path.insert(0, '.')
from fractions import Fraction
from quivermutation import arrowPaths as AP
def fullReduce(comb, pivots):
    row = {k: Fraction(v) for k, v in comb.items()}
    out = {}
    while row:
        head = min(row)
        if head in pivots:
            f = row[head]; pr = pivots[head]
            row = {k: row.get(k, Fraction(0)) - f * pr.get(k, Fraction(0)) for k in set(row) | set(pr)}
            row = {k: v for k, v in row.items() if v != 0}
        else:
            out[head] = row.pop(head)
    return out
if '--unpatched' in sys.argv: sys.argv.remove('--unpatched')
else: AP.reduceAgainstPivots = fullReduce
sys.argv = [sys.argv[0]] + sys.argv[1:]
import os
_src = open('workshop/rounds/014/scholar_walk.py').read()
_cap = int(os.environ.get('MAXEXP', '0'))   # MAXEXP=N: stop after N expansions (deterministic, replaces the time budget)
if _cap: _src = _src.replace('a.budget and time.time() - t0 > a.budget', 'expanded >= %d' % _cap)
exec(compile(_src, 'scholar_walk', 'exec'))
