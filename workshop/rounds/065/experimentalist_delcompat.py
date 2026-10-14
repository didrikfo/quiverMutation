"""S-1 item 1 at n = 9 -> 8: how often do same-class pairs stay same-class under removeVertex?
usage: python workshop/rounds/065/experimentalist_delcompat.py [N]   (N = 9 default)
Classes come from classification.classifyLength (n=9: ~4 min, n=8: ~40 s), cached in /tmp/s1del_A_<n>.csv."""
import sys, os, csv, collections
from quivermutation import classification as cl, nakayama, piecewiseHereditary as pwh, mutationClassTable as mct
N = int(sys.argv[1]) if len(sys.argv) > 1 else 9
def classes(n):
    f = '/tmp/s1del_A_%d.csv' % n
    if not os.path.exists(f):
        cl.classifyLength(n, fileName=f, printOutput=False, resume=True)
    t = mct.MutationClassTable.fromCSV(f, n)
    return {r[0]: r[1] for r in t.rows()}
hi, lo = classes(N), classes(N - 1)
print('n=%d: %d rows, %d classes; n=%d: %d rows, %d classes' % (N, len(hi), len(set(hi.values())), N-1, len(lo), len(set(lo.values()))))
img = {}   # (rowstring, i) -> class at N-1 or None
undef = collections.Counter()
for rs, c in hi.items():
    a = nakayama.LinearNakayamaAlgebra.fromRelationString(N, rs)
    for i in range(1, N + 1):
        r = pwh.removeVertex(N, list(a.relLengths), i)
        if r is None: img[rs, i] = None; undef[i] += 1; continue
        b = nakayama.LinearNakayamaAlgebra(r[0], r[1])
        img[rs, i] = lo[b.relationString()]
byclass = collections.defaultdict(list)
for rs, c in hi.items(): byclass[c].append(rs)
tot = sum(len(v)*(len(v)-1)//2 for v in byclass.values())
print('same-class pairs at n=%d: %d' % (N, tot))
print('undefined (non-admissible) images per i:', dict(sorted(undef.items())))
# same-vertex matrix: pairs (x,y) same class, deletion i on x, j on y (unordered pairs, x != y), both defined
def frac(i, j):
    same = both = 0
    for c, rows in byclass.items():
        cnt_i = collections.Counter(); cnt_j = collections.Counter()
        for rs in rows:
            a, b = img[rs, i], img[rs, j]
            if a is not None: cnt_i[a] += 1
            if b is not None: cnt_j[b] += 1
        ti, tj = sum(cnt_i.values()), sum(cnt_j.values())
        if i == j:
            both += ti*(ti-1)//2; same += sum(v*(v-1)//2 for v in cnt_i.values())
        else:  # ordered (x deleted at i, y deleted at j), x != y
            x_y = ti*tj - sum(1 for rs in rows if img[rs, i] is not None and img[rs, j] is not None)
            s = sum(cnt_i[k]*cnt_j[k] for k in cnt_i) - sum(1 for rs in rows if img[rs, i] is not None and img[rs, i] == img[rs, j])
            both += x_y; same += s
    return same, both
print('\nSAME vertex i on both: i, defined pairs, same-class pairs, fraction')
for i in range(1, N + 1):
    s, b = frac(i, i); print(i, b, s, '%.4f' % (s/b if b else float('nan')))
print('\nMATRIX fraction (x deleted at i, y at j), rows i, cols j')
for i in range(1, N + 1):
    print(i, ' '.join('%.3f' % (lambda t: t[0]/t[1] if t[1] else float('nan'))(frac(i, j)) for j in range(1, N + 1)))
# purity: for each class and i, number of image classes
print('\nper n=%d class: size; number of distinct n=%d image classes for each i (None-excluded)' % (N, N-1))
for c, rows in sorted(byclass.items(), key=lambda t: -len(t[1])):
    print('%-22s %4d ' % (c, len(rows)) + ' '.join(str(len({img[rs, i] for rs in rows} - {None})) for i in range(1, N + 1)))
# any single i with all classes pure?
pure = [i for i in range(1, N + 1) if all(len({img[rs, i] for rs in rows} - {None}) <= 1 for rows in byclass.values())]
print('vertices i with every class mapping into <= 1 image class:', pure)
# best uniform rule baseline: random same-class rate at n-1
import itertools
sizes = collections.Counter(lo.values()); T = len(lo)
print('baseline: P(two random distinct n=%d rows same class) = %.4f' % (N-1, sum(v*(v-1) for v in sizes.values())/(T*(T-1))))
# existence: for a same-class pair (x,y), is there SOME (i,j) with delta_i x ~ delta_j y ?  (S(x) = set of image classes)
print('\nEXISTS some (i,j): per class, pairs, pairs with S(x) & S(y) nonempty, pairs with i = j agreeing for some common i')
T = E = Ei = 0
for c, rows in sorted(byclass.items(), key=lambda t: -len(t[1])):
    S = [{img[rs, i] for i in range(1, N+1)} for rs in rows]; V = [[img[rs, i] for i in range(1, N+1)] for rs in rows]
    p = e = ei = 0
    for a in range(len(rows)):
        for b in range(a):
            p += 1; e += bool(S[a] & S[b]); ei += any(V[a][k] == V[b][k] for k in range(N))
    T += p; E += e; Ei += ei
    print('%-22s %7d %7d %7d' % (c, p, e, ei))
print('TOTAL', T, E, Ei, '%.4f %.4f' % (E/T, Ei/T))
