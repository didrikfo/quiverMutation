"""Round 135 revision (experimentalist): DIRECTED boundary census, answering the referee.  Repo root:
  timeout 10m .venv/bin/python workshop/rounds/135/experimentalist_boundary_dir.py N H
For each family mut and der, counts key-changing edges with exactly one end at K = 0, split by direction
(K=0 -> K>=1 and K>=1 -> K=0), with the head-shape tally (first 4 entries of source -> target) and whether the tail m[4:] is unchanged.
"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
import experimentalist_boundary as B
E = B.E

if __name__ == '__main__':
    n, h = int(sys.argv[1]), int(sys.argv[2])
    cls, _ = E.keyclass(n, h); ids = {}; lab = {m: ids.setdefault(E.key(n-1, m[1:]), len(ids)) for m in cls}
    K = {m: E.Kof(m, n) for m in cls}; s = lambda m: ''.join(map(str, m))
    for moves in ('mut', 'der'):
        fwd = collections.Counter(); rev = collections.Counter(); tailchg = 0; freeK = collections.Counter()
        for m in cls:
            for r, kd in B.kinds(n, m, moves):
                if r not in K or lab[r] == lab[m] or (K[m] == 0) == (K[r] == 0): continue
                shape = '%s->%s %s' % (s(m[:4]), s(r[:4]), kd)
                if K[m] == 0: fwd[shape] += 1; freeK[K[r]] += 1
                else: rev[shape] += 1
                tailchg += (m[4:] != r[4:])
        print('n=%d h=%d %s: K=0 -> K>=1: %d %s; free-end K %s | K>=1 -> K=0: %d %s | tail m[4:] changed on %d'
              % (n, h, moves, sum(fwd.values()), dict(fwd), dict(freeK), sum(rev.values()), dict(rev), tailchg))
