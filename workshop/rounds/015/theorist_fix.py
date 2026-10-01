"""Round 015 theorist: monkeypatch arrowPaths.reduceAgainstPivots to a FULL reduction (every term that is a pivot
column is eliminated, not just the leading one), then redo, for the replayed parents of E-084 n = 8 class 2
(1 rejection + 10 'M' lines), the mutation + Cartan congruence + key check.  Usage: theorist_fix.py [n class]"""
import sys
n, c = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (8, 2)
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
            out[head] = row.pop(head)      # leading term is not a pivot: keep it, continue with the tail
    return out
orig = AP.reduceAgainstPivots
if '--patched' in sys.argv or True:
    AP.reduceAgainstPivots = fullReduce
sys.argv = ['x', str(n), str(c)]
exec(open('workshop/rounds/015/theorist_cartan.py').read())
