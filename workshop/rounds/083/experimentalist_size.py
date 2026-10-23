"""S-1 class-level sizing (round 083). usage: experimentalist_size.py N   -> time, peak RSS, rows, orbits of freeMoves.derivedOrbits(N); saves orbit labels
to workshop/rounds/083/experimentalist_orbits_N.pkl only if --save is given (kept out of git)."""
import sys, time, resource, pickle, collections
from quivermutation import freeMoves
n = int(sys.argv[1]); t = time.time()
lnas, orbs = freeMoves.derivedOrbits(n)
dt = time.time() - t
sz = sorted((len(v) for v in orbs.values()), reverse=True)
print('n', n, 'rows', len(lnas), 'orbits', len(orbs), 'largest', sz[:6], 'singletons', sum(1 for s in sz if s == 1), '%.0fs' % dt, 'peakRSS %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if '--save' in sys.argv:
    lab = {}
    for i, v in enumerate(sorted(orbs.values(), key=lambda v: (-len(v), v[0]))):
        for r in v: lab[r] = i
    pickle.dump(lab, open('/tmp/claude-0/orb_%d.pkl' % n, 'wb'))
