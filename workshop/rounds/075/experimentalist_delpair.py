"""Round 075, S-1 item 2 on the E-190 pair family: A_n = kA_n/rad^(n-1) (row (n-1,0,..,0)) and B_n = one relation of
3 arrows on the last four vertices (row (0,..,0,3,0)).  For n = 5..N: Coxeter key of each, forward move orbit (free+edge+double
moves, mirror-closed) and whether A_n, B_n are joined by orbit + mirror; then every deletion of A_n and B_n: image row, image key
(equal to key(A_{n-1}) / key(B_{n-1}) ?), image orbit size and orbit-join status.
usage: timeout 10m .venv/bin/python workshop/rounds/075/experimentalist_delpair.py [Nmax=12]"""
import sys, time
from quivermutation import freeMoves as fm, coxeterTables as ct, piecewiseHereditary as pwh
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
def A(n): return tuple([n - 1] + [0] * (n - 3))
def B(n): return tuple([0] * (n - 4) + [3, 0]) if n >= 5 else None
def show(r): return ''.join(map(str, r))
def mir(r): return fm.mirrorRow(len(r) + 2, r)
def key(r): return ct.lnaCoxeterKey(len(r) + 2, r)
def kid(k, tab={}): return tab.setdefault(k, len(tab))
def orb(r, lim=200000):
    w = fm.orbitReport(len(r) + 2, r, free=True, limit=lim)
    return len(w.rows), w.stoppedBy, w.rows
def joined(a, b, lim=200000):
    n = len(a) + 2
    return (fm.movesJoin(n, a, b, free=True, limit=lim) is not None) or (fm.movesJoin(n, a, mir(b), free=True, limit=lim) is not None)
print('n | key(A)=key(B) | key id | orbit A | orbit B | joined by orbit+mirror')
for n in range(5, NMAX + 1):
    a, b = A(n), B(n); t = time.time()
    na, sa, _ = orb(a); nb, sb, _ = orb(b)
    print('%2d | %s | %d | %d %s | %d %s | %s  (%.1fs)' % (n, key(a) == key(b), kid(key(a)), na, sa, nb, sb, joined(a, b), time.time() - t))
    sys.stdout.flush()
print('\nDeletions: n, source, i, image row, image key == key(A_{n-1}) / key(B_{n-1}), image orbit size')
for n in range(6, NMAX + 1):
    for name, r in (('A', A(n)), ('B', B(n))):
        out = []
        for i in range(1, n + 1):
            d = pwh.removeVertex(n, list(r), i)
            img = tuple(d[1])
            out.append((i, show(img), key(img) == key(A(n - 1)), key(img) == key(B(n - 1)), orb(img, 50000)[0]))
        # compress runs of identical outcomes
        runs = []
        for o in out:
            if runs and tuple(runs[-1][1:]) == tuple(o[1:]): runs[-1][0].append(o[0])
            else: runs.append([[o[0]], *o[1:]])
        for idx, rimg, ka, kb, osz in runs:
            print('n=%2d %s i=%s -> %s  key=A_{n-1}:%s key=B_{n-1}:%s orbit %d' % (n, name, ('%d..%d' % (idx[0], idx[-1])) if len(idx) > 1 else idx[0], rimg, ka, kb, osz))
