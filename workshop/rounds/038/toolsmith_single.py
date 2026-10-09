"""Single-relation cores (length a, h free vertices at the head, K at the tail, h+K+a+1=n): Coxeter key vs offset.
E-127 / the 037 mirror table says the I1/I2 split at n=11 is carried by a single 3 (or 7) core at K_eff=3 vs 4.
Round 038 toolsmith copy of rounds/037/maverick_single.py with the label block fixed (label dicts are keyed by rows of length n-2, not n).
usage: python workshop/rounds/038/toolsmith_single.py"""
import collections
from quivermutation import coxeterTables as ct
def key(n, a, h):
    row = [0] * n; row[h] = a      # relation starting at vertex h+1 (0-based position h), length a arrows
    return ct.lnaCoxeterKey(n, tuple(row))
for a in (3, 4, 7):
    for n in (9, 10, 11, 12):
        ks = collections.defaultdict(list)
        for h in range(0, n - a):
            K = n - a - 1 - h
            ks[key(n, a, h)].append((h, K))
        print('a', a, 'n', n, 'classes by key:', [sorted(v) for v in ks.values()])

print('--- minimal failing n for threshold K0 (single 3: h != K, both >= K0), by key at n-1')
for K0 in (3, 4, 5):
    for n in range(8, 20):
        a = 3
        bad = []
        for h in range(K0, n - a - 1 - K0 + 1):
            K = n - a - 1 - h
            if K < K0 or h >= K: continue
            if key(n - 1, a, h) != key(n - 1, a, h - 1 + 0) and key(n - 1, a, h) is not None:
                pass
            t1 = key(n - 1, a, h)          # tail deleted from (h,K): (h,K-1)
            t2 = key(n - 1, a, h - 1)      # head deleted: (h-1,K)
            if t1 != t2: bad.append((h, K))
        if bad: print('K0', K0, 'first n with a failing single-3 class:', n, bad); break
print('--- identify I1/I2 at n=10 via the 034 labels (keys are rows of length n-2)')
import sys
sys.path.insert(0, "workshop/rounds/030"); sys.path.insert(0, "workshop/rounds/027")
import maverick_endstrip2 as e2
def labelOfRow(lab, n, row):
    """label of a row given at length n (relation start index h, length a); label dicts use length n-2"""
    r = tuple(row)
    if r not in lab: r = r[:n - 2]
    assert all(x == 0 for x in row[n - 2:]) and r in lab, ('row not a valid LNA', row)
    return lab[r]
if __name__ == '__main__':
    lab2, _ = e2.corrected(10)
    for h in (4, 3, 2):
        row = [0] * 10; row[h] = 3
        print('n=10 single 3 (h,K)=', (h, 10 - 4 - h), labelOfRow(lab2, 10, row))
