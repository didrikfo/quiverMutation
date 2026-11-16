"""Round 129 (maverick): does the LIBRARY walk itself (`search.mutationSearchDepthFirst`, F-only, coxeterGuard=True, default gate,
fingerprint.Visited dedup) reach a node whose SNF(C + C^T) differs from the root's?  Independent of the round-127/129 BFS code:
only the library walk and a visitor that computes SNF(S) (python-flint) and the signature of S.
usage: [REFUSE=1] timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_libwalk.py N DEPTH ROOT [dual]
(REFUSE=1 passes refuseNonTilting=True: the control on which the crossing must disappear)"""
import sys, os, time, collections
sys.path.insert(0, '.')
import numpy as np, flint
from quivermutation import search, fingerprint, invariants, pathAlgebra, nakayama as nk
N, D, rt = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
side = sys.argv[4] if len(sys.argv) > 4 else 'root'
def inv(A):
    C = invariants.integerCartanMatrix(A); n = len(C)
    S = [[int(C[i][j]) + int(C[j][i]) for j in range(n)] for i in range(n)]
    M = flint.fmpz_mat(S).snf(); d = tuple(sorted(abs(int(M[i, i])) for i in range(n)))
    ev = np.linalg.eigvalsh(np.array(S, dtype=float)); return d, (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
r = nk.LinearNakayamaAlgebra.fromClassName(rt); start = r if side == 'root' else pathAlgebra.dualPathAlgebra(r)
i0 = inv(start); T = collections.Counter(); hits = []; t0 = time.time()
def visitor(A, path):
    T['visits'] += 1
    i = inv(A)
    if i != i0:
        hits.append((tuple(path), i[0][-3:], i[1]))
        if len(hits) <= 5: print('  OTHER SNF CLASS: path %s SNF tail %s sig %s' % (' '.join(map(str, path)), i[0][-3:], i[1]), flush=True)
search.mutationSearchDepthFirst(start, D, printOutput=False, visitor=visitor, visited=fingerprint.Visited(), refuseNonTilting=os.environ.get('REFUSE') == '1')
print('LIBWALK%s %s %s n=%d depth %d: root SNF tail %s sig %s; visits %d; nodes in another (SNF, sig) class %d  %.0fs'
      % (' refuseNonTilting' if os.environ.get('REFUSE') == '1' else '', rt, side, N, D, i0[0][-3:], i0[1], T['visits'], len(hits), time.time() - t0))
