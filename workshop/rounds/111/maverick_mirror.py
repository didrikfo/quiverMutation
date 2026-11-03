"""Round 111 response: basis for dropping tail ends (mirror equivariance) and E-173 reconcile.
usage: maverick_mirror.py                # all n=11 rows (exhaustive), 600 random n=13 rows, n=17 K1=6 class (head rows of the lone-3 key class)
Tests: (a) lnaCoxeterKey(r)==lnaCoxeterKey(mirror r); (b) key(removeVertex(n,r,n)) == key(removeVertex(n,mirror r,1)); None agree.
(c) at n=17: for every head-end class row m (K>=6), image key of tail deletion of mirror(m) equals image key of head deletion of m."""
import random, sys, time
from quivermutation import coxeterTables as ct, nakayama, piecewiseHereditary as pwh, freeMoves as fm
def lone(n, h): return tuple([0]*h + [3] + [0]*(n-2-h-1))
def img(n, r, v):
    x = pwh.removeVertex(n, list(r), v)
    return None if x is None else ct.lnaCoxeterKey(n-1, tuple(x[1]))
def test(n, rows, tag):
    bad_key = bad_img = bad_none = nn = 0
    for r in rows:
        r = tuple(r); m = fm.mirrorRow(n, r)
        if ct.lnaCoxeterKey(n, r) != ct.lnaCoxeterKey(n, m): bad_key += 1
        a, b = img(n, r, n), img(n, m, 1)
        if (a is None) != (b is None): bad_none += 1
        elif a is None: nn += 1
        elif a != b: bad_img += 1
    print(tag, 'rows', len(rows), 'key-mirror failures', bad_key, 'image-key failures', bad_img, 'None mismatches', bad_none, 'both None', nn, flush=True)
t0 = time.time()
A = [tuple(r) for r in nakayama.allRelationLengths(11)]
test(11, A, 'n=11 exhaustive')
random.seed(1)
B = [tuple(r) for r in nakayama.allRelationLengths(13)]
test(13, random.sample(B, 600), 'n=13 random')
n, K1 = 17, 6
h0 = (n-4)//2
tg = {ct.lnaCoxeterKey(n, lone(n, h)) for h in (h0, n-4-h0)}
rows = {(0,)*K1+tuple(s) for s in nakayama.allRelationLengths(n-K1)}
cls = [r for r in rows if ct.lnaCoxeterKey(n, r) in tg]
test(n, cls, 'n=17 K1=6 key class head rows')
# E-173 reconcile: lone rows (6,7) and (7,6) at n=17; head vs tail deletion
for h in (6, 7):
    r = lone(n, h); x1 = pwh.removeVertex(n, list(r), 1); xn = pwh.removeVertex(n, list(r), n)
    print('lone h=%d: head image row %s | tail image row %s | row equal %s | key equal %s' % (h, tuple(x1[1]), tuple(xn[1]), tuple(x1[1]) == tuple(xn[1]),
          ct.lnaCoxeterKey(16, tuple(x1[1])) == ct.lnaCoxeterKey(16, tuple(xn[1]))))
print('%.0fs' % (time.time()-t0))
