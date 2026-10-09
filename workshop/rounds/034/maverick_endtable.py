"""Per-end table of the n=11 K>=3 free-end failure (E-120). For the failing class (1305 LNAs): every deletable end
(free run K>=3), its side, K, core length L = last-first+1 (+ gap pattern), source orbit, image class.
usage: python workshop/rounds/034/maverick_endtable.py"""
import sys, collections
sys.path.insert(0, "workshop/rounds/030"); sys.path.insert(0, "workshop/rounds/027")
import maverick_endstrip2 as e2
from quivermutation import piecewiseHereditary as pwh, freeMoves
n = 11; lab, _ = e2.corrected(n); lab2, _ = e2.corrected(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
cls = collections.defaultdict(list)
for m, l in lab.items(): cls[l].append(m)
tgt = [c for c in cls if str(c).endswith('-3, -3, -2, -1, 0, 1, 1)') and len(cls[c]) == 1305][0]
imgs = {}; 
def name(t):
    k = imgs.setdefault(t, 'I%d' % (len(imgs) + 1)); return k
lnas, orbits = freeMoves.derivedOrbits(n, rules=None, free=True, edges=True, doubles=True)
oid = {m: o for o, mem in orbits.items() for m in mem}
rows = []
for m in cls[tgt]:
    R = sp(m); first = min(s for s, e in R); last = max(e for s, e in R)
    for side, K, v in (('head', first - 1, 1), ('tail', n - last, n)):
        if K >= 3:
            t = lab2[tuple(pwh.removeVertex(n, list(m), v)[1])]
            rows.append((side, K, oid[m], name(t), last - first + 1, ''.join(map(str, m)), sorted(R)))
print('class', tgt); print('images', {v: str(k) for k, v in imgs.items()})
c = collections.Counter((r[1], r[0], r[2], r[3]) for r in rows)
print('K side orbit image count')
for k in sorted(c): print(*k, c[k])
print('--- by (K, core length L, image)')
c = collections.Counter((r[1], r[4], r[3]) for r in rows)
for k in sorted(c): print(*k, c[k])
print('--- by (K, free run at OTHER end, image)')
c = collections.Counter()
for r in rows:
    R = r[6]; first = min(s for s, e in R); last = max(e for s, e in R)
    other = (n - last) if r[0] == 'head' else first - 1
    c[(r[1], other, r[3])] += 1
for k in sorted(c): print(*k, c[k])
print('--- core words (a_i string) with K=3 per image, 6 examples each')
for I in sorted(set(r[3] for r in rows)):
    ex = [(r[0], r[5]) for r in rows if r[3] == I and r[1] == 3][:6]; print(I, ex)
print('distinct (core word stripped of leading/trailing zeros) by image for K=3:')
d = collections.defaultdict(set)
for r in rows:
    if r[1] == 3: d[r[3]].add(r[5].strip('0') if r[0] else r[5])
for I in d: print(I, len(d[I]), sorted(d[I])[:8])
print('=== is the image a function of (K, stripped core word oriented away from the deleted end)?')
f = collections.defaultdict(set)
for r in rows:
    w = r[5]; w = w if r[0] == 'tail' else w[::-1]   # word read so that the deleted end is on the right (arrow-index reversal is a guess; mirrorRow would be exact)
    f[(r[1], w.strip('0'))].add(r[3])
print('keys', len(f), 'with 2 images', sum(len(v) > 1 for v in f.values()))
print('K=3 I1 only:', sorted(k[1] for k, v in f.items() if k[0] == 3 and v == {'I1'}))
print('K=3 I2 only:', sorted(k[1] for k, v in f.items() if k[0] == 3 and v == {'I2'}))
# image vs. the free zero-padding inside the stripped part (interior zeros / number of zeros at far end)
print('=== same oriented core word at K=3 and K=4?')
w3 = {k[1]: v for k, v in f.items() if k[0] == 3}; w4 = {k[1]: v for k, v in f.items() if k[0] == 4}
print('K=4 words', sorted(w4)); print('common', {w: (w3[w], w4[w]) for w in w3 if w in w4})
# per-orbit count of K=3 ends and distinct cores
print('orbit sizes in class', collections.Counter(oid[m] for m in cls[tgt]).most_common(5), 'orbits', len(set(oid[m] for m in cls[tgt])))
