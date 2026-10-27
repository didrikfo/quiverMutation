"""Round 091 (experimentalist), T10 calibration: scan EVERY node of a ball (state pickle of toolsmith_d6scan.py build) for key-kept perI (J != 0) steps.
usage (repository root): experimentalist_jscan.py STATE.pkl OUTDIR [MAXDEPTH]
Writes OUTDIR/j<k>.pkl in the skeptic_kept.py format (one dict per file) and prints one line per kept step: node depth, path, kind, vertex."""
import sys, pickle, os
A_ = sys.argv; ST, OUT = A_[1], A_[2]; MAXD = int(A_[3]) if len(A_) > 3 else 99
sys.argv = ['x', 'ctl', 'x']
src = open('workshop/rounds/087/toolsmith_d6scan.py').read().split("if MODE == 'build'")[0]
sys.path.insert(0, '.'); exec(compile(src, 's', 'exec'))
S = pickle.load(open(ST, 'rb')); seen = S['seen']; os.makedirs(OUT, exist_ok=True)
def path(k):
    p = []
    while seen[k][1] is not None: p.append(seen[k][2]); k = seen[k][1]
    return p[::-1]
import collections; cnt = collections.Counter(); n = 0
for k, (x, pk, mv, d) in seen.items():
    if d > MAXD: continue
    cnt['nodes_d%d' % d] += 1
    for kind, v, st in jcands(x):
        cnt['J:' + st] += 1
        if st == 'kept':
            p = path(k)
            pickle.dump([dict(parent=x, kind=kind, v=v, depth=d, path=p, rootalg=None)], open('%s/j%d.pkl' % (OUT, n), 'wb'))
            print('KEPT j%d depth %d path %s %s %d' % (n, d, p, kind, v), flush=True); n += 1
print(dict(cnt))
