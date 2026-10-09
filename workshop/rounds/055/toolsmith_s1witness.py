"""S-1 n = 15, K0 = 5: witness check without the full class scan (round 055).
usage: python workshop/rounds/055/toolsmith_s1witness.py
Lone 3 at (h,K) has h+K = n-4. Uses piecewiseHereditary.removeVertex (as E-144 'ends') and lnaCoxeterKey."""
import sys
from quivermutation import coxeterTables as ct, piecewiseHereditary as pwh, nakayama
def lone(n, h): return tuple([0]*h + [3] + [0]*(n-2-h-1))
def key(n, row): return ct.lnaCoxeterKey(n, tuple(row))
for n, hs in ((11, (3, 4)), (13, (4, 5)), (15, (5, 6)), (15, (4, 7)), (16, (6, 6)), (17, (6, 7))):
    h1, h2 = hs
    a, b = lone(n, h1), lone(n, h2)
    print('n', n, 'lone3 (h,K)', (h1, n-4-h1), 'vs', (h2, n-4-h2), 'same key:', key(n, a) == key(n, b))
    for h in hs:
        m = lone(n, h); K = n-4-h
        head = tuple(pwh.removeVertex(n, list(m), 1)[1]); tail = tuple(pwh.removeVertex(n, list(m), n)[1])
        print('   (h,K)=(%d,%d) head-del %s tail-del %s keys equal: %s' % (h, K, head.index(3), tail.index(3), key(n-1, head) == key(n-1, tail)))
# direct: all single-3 images at n=15, key at n-1
n = 15
cls = {}
for h in range(0, n-3):
    cls.setdefault(key(n, lone(n, h)), []).append((h, n-4-h))
print('n = 15 single-3 key classes:', list(cls.values()))
