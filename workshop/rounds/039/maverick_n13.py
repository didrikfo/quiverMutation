"""S-1 at n = 13: lone 3 at (h,K) = (4,5), (5,4). Image keys of every free-end deletion, then a bounded move orbit.
usage: python workshop/rounds/039/maverick_n13.py [orbitLimit]   (orbitLimit 0 = skip orbit)"""
import sys, time
from quivermutation import coxeterTables as ct, freeMoves, piecewiseHereditary as pwh

def row(n, a, h):
    r = [0] * n; r[h] = a; return tuple(r)
def key(n, a, h): return ct.lnaCoxeterKey(n, row(n, a, h))
def rem(n, r, v):
    return tuple(pwh.removeVertex(n, list(r), v)[1])

n = 13
for (h, K) in ((4, 5), (5, 4)):
    r = row(n, 3, h)
    print('source (h,K)=', (h, K), 'row', ''.join(map(str, r)), 'key', key(n, 3, h))
    for side, v, ok in (('head', 1, h), ('tail', n, K)):
        im = rem(n, r, v)
        print('  delete', side, 'free', ok, '-> n=12 row', ''.join(map(str, im)), 'key', ct.lnaCoxeterKey(n - 1, im))
print('keys n=13 (4,5) vs (5,4) equal:', key(13, 3, 4) == key(13, 3, 5))
print('n=12 images: (4,4)', key(12, 3, 4), ' (3,5)', key(12, 3, 3), ' equal:', key(12, 3, 4) == key(12, 3, 3))
lim = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
if lim:
    for h in (4, 5):
        t = time.time()
        rep = freeMoves.orbitReport(n, list(row(n, 3, h)[:n - 2]), free=True, edges=True, doubles=True, limit=lim,
                                    target=row(n, 3, 9 - h)[:n - 2])
        print('orbit from h =', h, 'stoppedBy', rep.stoppedBy, 'rows', len(rep), 'found', rep.found, '%.0fs' % (time.time() - t))
