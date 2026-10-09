"""Round 045: guarded depth-D BFS from every LNA of length n (no dual); tally gate+key-admitted edges by tiltingPlus (J = 0 or not).
This is the walk shape behind mergeReport / E-031 depth-4 merges.  Usage: .venv/bin/python workshop/rounds/045/experimentalist_t10b_sweep.py N DEPTH [SECS]"""
import sys, time
sys.argv, args = ['x', '--plan'], sys.argv[1:]
src = open('workshop/rounds/045/experimentalist_t10b.py').read().split("P = pairs()")[0]
exec(compile(src, 't10b', 'exec'))
n, D = int(args[0]), int(args[1]); secs = float(args[2]) if len(args) > 2 else 500
rows = [r for r in __import__('itertools').product(*[range(0, 1)]*0)]
import batch
from quivermutation import lnaMoves
lnas = [tuple(r) for r in (fm.derivedOrbits(n, rules=lm.ALL_MOVES, free=False, edges=True, doubles=True)[0])]
print('n', n, 'LNAs', len(lnas), 'depth', D, flush=True)
tot = Counter(); T0 = time.time(); done = 0; bad = []
for r in lnas:
    if time.time() - T0 > secs: break
    st = Counter()
    reach(nk.LinearNakayamaAlgebra(n, list(r)), D, 'guard', st)
    tot.update(st); done += 1
    if st['J!=0']: bad.append((''.join(map(str, r)), st['J!=0']))
print('done %d of %d LNAs in %.0fs (%s); edges %s; LNAs with J!=0 admitted edges: %s' % (done, len(lnas), time.time()-T0, 'complete' if done == len(lnas) else 'CAPPED', dict(tot), bad), flush=True)
