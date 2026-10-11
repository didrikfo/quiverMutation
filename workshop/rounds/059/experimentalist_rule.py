"""Test the predicate (see pred): right mutation at v of an LNA r (r_j = arrows of the relation starting at j, 0 = none) is an LNA.
Reads experimentalist_which_<n>.json written by experimentalist_which.py.  Usage: rule.py N"""
import sys, json, collections
n = int(sys.argv[1])
rows = json.load(open(f"workshop/rounds/059/experimentalist_which_{n}.json"))
def pred(r, v):
    # relation (j, m) covers arrows j..j+m-1.  R(v) iff arrow v is covered by no relation starting before v,
    # and (r_v == 2, or v == n-1 where a virtual length-2 relation would stick out of the line)
    if v > n - 1: return False
    for j in range(1, v):
        if j <= n - 2 and r[j - 1] and j + r[j - 1] - 1 >= v: return False
    return v == n - 1 or r[v - 1] == 2
bad = collections.Counter(); tot = 0
for r, res in rows:
    actual = {v for v, _ in res.get('R', [])}
    pr = {v for v in range(1, n + 1) if pred(r, v)}
    tot += len(actual)
    for v in actual ^ pr: bad[(v, v in actual)] += 1
print("n", n, "R events", tot, "mismatches", dict(bad))
# image check: new vector = r with (r_{v-1}, r_v) := (2, 0)  (v = 1: r_1 := 0 only; v = n-1: r_{n-2} := 2 only)
def image(r, v):
    r = list(r)
    if v <= n - 2: r[v - 1] = 0
    if v >= 2: r[v - 2] = 2
    return r
badimg = 0; seen = 0
for r, res in rows:
    for v, new in res.get('R', []):
        seen += 1
        if image(r, v) != new: badimg += 1
print("image mismatches", badimg, "of", seen)
# left side by mirror: opposite algebra reverses the line; test L(v) <-> R(n-v) on the mirrored vector
def mirror(r):
    out = [0] * (n - 2)
    for j in range(1, n - 1):
        m = r[j - 1]
        if m:
            out[n - j - m + 1 - 1] = m   # arrows j..j+m-1 -> n-(j+m-1)..n-j ; start n-j-m+1
    return out
for c in range(n - 1, n + 4):      # try L(v) <-> R(c - v) on the mirrored vector
    bad = 0
    for r, res in rows:
        actual = {v for v, _ in res.get('L', [])}
        mr = mirror(r)
        pr = {c - w for w in range(1, n) if pred(mr, w)}
        bad += len(actual ^ pr)
    print("left offset c =", c, "mismatching (LNA, v) pairs:", bad)
