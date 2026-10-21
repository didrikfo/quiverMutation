"""Round 078 response: (ctrl) positive control for experimentalist_join10.py; (targets) how many LNAs/duals of the target key lie in the balls.
  ctrl    random F/R walks of length L from 90000000, then join walk-end X back to the root with the same machinery (known path <= L).
  targets canonicalKey of every n=10 LNA and dual with the Coxeter key of 90000000, intersected with the A-ball (5) and B-ball (4).
Run from repo root; imports the machinery by exec'ing the top of experimentalist_join10.py."""
import sys, random, time
MY = sys.argv
src = open('workshop/rounds/078/experimentalist_join10.py').read().split("if MYARGV[1] == 'plan':")[0]
sys.argv = ['x', 'run']; exec(compile(src, 'j10', 'exec'))
if MY[1] == 'ctrl':
    random.seed(int(MY[2]))
    for L, da, dx in ((6, 3, 3), (8, 4, 4)):
        for trial in range(2 if L == 8 else 3):
            x = A; path = [x]
            for _ in range(L):
                ms = moves(x); kind, v, x = random.choice(ms); path.append(x)
            ends = fingerprint.canonicalKey(x)
            print('walk L=%d trial %d end key == root key: %s' % (L, trial, ends == fingerprint.canonicalKey(A)), flush=True)
            sa, _ = sideBall(A, da, 100000, 'A'); sx, _ = sideBall(x, dx, 100000, 'X')
            common = set(sa) & set(sx)
            print('CTRL L=%d (%d+%d) trial %d: shared %d' % (L, da, dx, trial, len(common)),
                  'HIT depth', min(sa[k] + sx[k] for k in common) if common else 'MISS', flush=True)
else:
    from quivermutation import pathAlgebra
    t = time.time(); T = []
    for r in nk.allRelationLengths(N):
        if ct.lnaCoxeterKey(N, tuple(r)) == kA:
            a = alg(r); T.append((tuple(r), 'lna', fingerprint.canonicalKey(a))); T.append((tuple(r), 'dual', fingerprint.canonicalKey(pathAlgebra.dualPathAlgebra(a))))
    print('targets', len(T), 'distinct canonicalKeys', len({k for _, _, k in T}), 'None keys', sum(k is None for _, _, k in T), '%.0fs' % (time.time() - t), flush=True)
    sa, _ = sideBall(A, 5, 30000, 'A'); sb, _ = sideBall(B, 4, 30000, 'B')
    for nm, s in (('A-ball r5', sa), ('B-ball r4', sb)):
        hit = [(r, w) for r, w, k in T if k in s]
        print(nm, 'target LNAs/duals inside:', len(hit), hit[:10])
    for nm, root in (('root A', A), ('root B', B)):
        print(nm, 'is target:', any(k == fingerprint.canonicalKey(root) for _, _, k in T))
