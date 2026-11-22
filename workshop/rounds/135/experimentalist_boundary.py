"""Round 135 (experimentalist): the key-changing BOUNDARY edges (one end K = 0, other K >= 1) of experimentalist_bridge.py.  Repo root:
  timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary.py N H [moves]   (default mut)
Prints, per K of the free end, the count of key-changing boundary edges, the move kind (rules / edge / double / strip), and the
distinct (K=0 row, K>=1 row) pairs (up to 12), with the free end's image key.
"""
import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '123'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '117'))
import experimentalist_classlevel as E
from quivermutation import lnaMoves, edgeMoves, doubleMutation, freeMoves


def kinds(n, m, moves):
    out = [(tuple(x), 'rules') for x in lnaMoves.rewritesOf(n, list(m))]
    if moves in ('mut', 'der'):
        out += [(tuple(x), 'edge') for x, _ in edgeMoves.rewritesOf(n, m)]
        out += [(tuple(x), 'double') for x, _ in doubleMutation.rewritesOf(n, m)]
    if moves == 'der': out.append((tuple(freeMoves.stripLengthTwo(m)), 'strip'))
    return out


if __name__ == '__main__':
    n, h = int(sys.argv[1]), int(sys.argv[2]); moves = sys.argv[3] if len(sys.argv) > 3 else 'mut'
    cls, _ = E.keyclass(n, h); ids = {}; lab = {m: ids.setdefault(E.key(n-1, m[1:]), len(ids)) for m in cls}
    K = {m: E.Kof(m, n) for m in cls}; byK = collections.Counter(); kind = collections.Counter(); pairs = set(); zfirst = collections.Counter()
    for m in cls:
        if K[m] != 0: continue
        for r, kd in kinds(n, m, moves):
            if K[r] >= 1 and lab[r] != lab[m]:
                byK[K[r]] += 1; kind[kd] += 1; pairs.add((m, r)); zfirst[m[0]] += 1
    s = lambda m: ''.join(map(str, m))
    print('n=%d h=%d K0=%d moves=%s: key-changing boundary edges (from the K=0 side) %d; by K of free end %s; by kind %s; by first entry of K=0 row %s'
          % (n, h, min(h, n-4-h), moves, sum(byK.values()), dict(sorted(byK.items())), dict(kind), dict(zfirst)))
    for m, r in sorted(pairs)[:12]:
        print('   %s [img %d] -> %s [K%d img %d]' % (s(m), lab[m], s(r), K[r], lab[r]))
