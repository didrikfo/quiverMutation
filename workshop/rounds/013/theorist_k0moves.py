"""One-step reduced moves out of k@0 (and the size of its orbit), labelled.  .venv/bin/python workshop/rounds/013/theorist_k0moves.py N kmin kmax"""
import sys
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/013')
import batch
from quivermutation import freeMoves as fm
from theorist_path import labelled
n = int(sys.argv[1])
for k in range(int(sys.argv[2]), int(sys.argv[3]) + 1):
    r = batch._rowFor(n, str(k), 0)
    rep = fm.orbitReport(n, tuple(r), free=fm.REDUCED, limit=300000)
    print('n=%d %d@0 orbit %d (%s)' % (n, k, len(rep.rows), rep.stoppedBy))
    for row, lab in labelled(n, tuple(r)).items(): print('   ', ''.join(map(str, row)), '<-', lab)
