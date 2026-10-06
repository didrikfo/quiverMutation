"""E-125 redone with freeMoves.mirrorRow for head ends (r034 reversed digit strings). Reads every K>=3 end of the failing
n=11 class as a TAIL word: tail ends keep m, head ends use mirrorRow(n, m). Reports whether image = f(K, tail word),
head-vs-tail agreement, and the words common to K=3 and K=4. usage: python workshop/rounds/037/maverick_endtable_mirror.py"""
import sys, collections
sys.path.insert(0, "workshop/rounds/030"); sys.path.insert(0, "workshop/rounds/027")
import maverick_endstrip2 as e2
from quivermutation import piecewiseHereditary as pwh, freeMoves
n = 11; lab, _ = e2.corrected(n); lab2, _ = e2.corrected(n - 1)
def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
cls = collections.defaultdict(list)
for m, l in lab.items(): cls[l].append(m)
tgt = [c for c in cls if str(c).endswith('-3, -3, -2, -1, 0, 1, 1)') and len(cls[c]) == 1305][0]
imgs = {}
def name(t): return imgs.setdefault(t, 'I%d' % (len(imgs) + 1))
rows = []
for m in cls[tgt]:
    R = sp(m); first = min(s for s, e in R); last = max(e for s, e in R)
    for side, K, v in (('head', first - 1, 1), ('tail', n - last, n)):
        if K >= 3:
            t = lab2[tuple(pwh.removeVertex(n, list(m), v)[1])]
            w = tuple(m) if side == 'tail' else freeMoves.mirrorRow(n, m)
            rows.append((side, K, name(t), w))
print('images', {v: str(k) for k, v in imgs.items()}); print('ends', len(rows))
def word(w):  # tail-oriented: core only, as (first-1 zeros stripped) digit string
    return ''.join(map(str, w)).strip('0')
f = collections.defaultdict(set); fs = collections.defaultdict(set)
for side, K, I, w in rows:
    f[(K, word(w))].add(I); fs[(K, word(w), side)].add(I)
print('keys', len(f), 'with 2 images', sum(len(v) > 1 for v in f.values()))
print('head/tail disagree on same (K,word):', sum(1 for k in f if ('head' in {s for (kk, ww, s) in fs if (kk, ww) == k}) and ('tail' in {s for (kk, ww, s) in fs if (kk, ww) == k}) and len({i for (kk, ww, s), v in fs.items() if (kk, ww) == k for i in v}) > 1))
print('head-only keys', sum(1 for k in f if not any((k[0], k[1], 'tail') == x for x in fs)), 'tail-only', sum(1 for k in f if not any((k[0], k[1], 'head') == x for x in fs)))
for K in (3, 4):
    for I in sorted(set(i for (kk, w), v in f.items() if kk == K for i in v)):
        print('K', K, I, sorted(w for (kk, w), v in f.items() if kk == K and v == {I}))
w3 = {k[1]: v for k, v in f.items() if k[0] == 3}; w4 = {k[1]: v for k, v in f.items() if k[0] == 4}
print('common', {w: (sorted(w3[w]), sorted(w4[w])) for w in w3 if w in w4})
print('K=3 words ending in 3 -> images', collections.Counter(tuple(sorted(v)) for w, v in w3.items() if w.endswith('3')))
print('K=3 words not ending in 3 -> images', collections.Counter(tuple(sorted(v)) for w, v in w3.items() if not w.endswith('3')))

print('=== free-move normal form: length-2 relations are free, so drop them; K_eff = free run beyond the last NON-2 relation')
g = collections.defaultdict(set); ex = {}
cnt = collections.Counter()
for m in cls[tgt]:
    for side, v in (('head', 1), ('tail', n)):
        w = tuple(m) if side == 'tail' else freeMoves.mirrorRow(n, m)
        R = [(i + 1, i + 1 + a) for i, a in enumerate(w) if a]
        first = min(s for s, e in R); last = max(e for s, e in R)
        if n - last < 3: continue
        t = lab2[tuple(pwh.removeVertex(n, list(m), v)[1])]
        N2 = [(s, e) for s, e in R if e - s != 2]
        if not N2: key = ('only2', n - last); 
        else:
            lastN = max(e for s, e in N2); f0 = min(s for s, e in N2)
            key = (n - lastN, tuple((s - f0, e - s) for s, e in sorted(N2)))
        g[key].add(name(t)); cnt[(key, name(t))] += 1
print('keys', len(g), 'with 2 images', sum(len(v) > 1 for v in g.values()))
for k in sorted(g, key=lambda k: (str(k[0]), len(k[1]) if len(k) > 1 else 0, str(k))):
    print(k, sorted(g[k]))
