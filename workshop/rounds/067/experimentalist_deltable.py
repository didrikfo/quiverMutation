"""S-1 item 2 at n = 9 -> 8: table of image classes under deletion of vertex i, and which invariant changes by a rule.
usage: .venv/bin/python workshop/rounds/067/experimentalist_deltable.py   (reads the round-065 cached class CSVs; ~13 s)"""
import collections, csv, sys
from quivermutation import nakayama, piecewiseHereditary as pwh, coxeterTables as ct
N = 9
D = 'workshop/rounds/065/'
def load(n):
    return {r['relations']: r['class'] for r in csv.DictReader(open(D + 'experimentalist_classes_n%d.csv' % n))}
hi, lo = load(N), load(N - 1)
rows = {}
for rs, c in hi.items():
    a = nakayama.LinearNakayamaAlgebra.fromRelationString(N, rs); rows[rs] = (c, tuple(a.relLengths))
cls = sorted({c for c, _ in rows.values()}, key=lambda c: (-sum(1 for v in rows.values() if v[0] == c), c))
lcls = sorted(set(lo.values()), key=lambda c: (-sum(1 for v in lo.values() if v == c), c))
L = {c: 'abcdefghijk'[k] for k, c in enumerate(lcls)}
print('image classes (n=8):', {L[c]: (c, sum(1 for v in lo.values() if v == c)) for c in lcls})
img = {}; key = {}
for rs, (c, rl) in rows.items():
    key[rs] = ct.lnaCoxeterKey(N, rl)
    for i in range(1, N + 1):
        r = pwh.removeVertex(N, list(rl), i); b = nakayama.LinearNakayamaAlgebra(r[0], r[1])
        img[rs, i] = (lo[b.relationString()], ct.lnaCoxeterKey(N - 1, tuple(r[1])))
by = collections.defaultdict(list)
for rs, (c, _) in rows.items(): by[c].append(rs)
print('\nTABLE 1: source class (size); per vertex i=1..9: image-class histogram a:count,... (letters = n=8 classes by size)')
for c in cls:
    print('%s (%d)' % (c, len(by[c])))
    for i in range(1, N + 1):
        h = collections.Counter(L[img[rs, i][0]] for rs in by[c])
        print('   i=%d  %s' % (i, ' '.join('%s:%d' % kv for kv in sorted(h.items(), key=lambda t: -t[1]))))
# key facts
print('\nTABLE 2: per source class, distinct source keys; per i distinct image classes / image keys')
for c in cls:
    ks = {key[rs] for rs in by[c]}
    print('%-22s %4d srckeys=%d  ' % (c, len(by[c]), len(ks)) + ' '.join('%d/%d' % (len({img[rs, i][0] for rs in by[c]}), len({img[rs, i][1] for rs in by[c]})) for i in range(1, N + 1)))
# local features of vertex i in a row
def feat(rl, i):
    # relation starting at v with r arrows covers arrows v..v+r-1 (arrow j: vertex j -> j+1)
    cov = lambda arrow: sum(1 for v, r in enumerate(rl, 1) if r and v <= arrow < v + r)
    cl = cov(i - 1) if i > 1 else None; cr = cov(i) if i < N else None
    start = rl[i - 1] if i <= len(rl) else 0
    end = sum(1 for v, r in enumerate(rl, 1) if r and v + r == i)
    both = sum(1 for v, r in enumerate(rl, 1) if r and v <= i - 1 and i <= v + r - 1)
    return ('cl%s' % cl, 'cr%s' % cr, 'both%d' % both)
def determinism(name, f, target):
    """fraction of classes (and of (row,i) cells) for which target is a function of (source class, f(row,i))"""
    g = collections.defaultdict(collections.Counter)
    for rs, (c, rl) in rows.items():
        for i in range(1, N + 1):
            k = f(c, rl, i)
            if k is not None: g[k][target(rs, i)] += 1
    tot = sum(sum(h.values()) for h in g.values()); ok = sum(max(h.values()) for h in g.values())
    print('%-46s groups=%4d  cells=%6d  modal share=%.4f  pure groups=%d' % (name, len(g), tot, ok / tot, sum(1 for h in g.values() if len(h) == 1)))
tc = lambda rs, i: img[rs, i][0]; tk = lambda rs, i: img[rs, i][1]
print('\nTABLE 3: is the image a function of (source class, local feature of i)?')
R = {rs: rl for rs, (c, rl) in rows.items()}; C = {rs: c for rs, (c, rl) in rows.items()}
for nm, f in [('class', lambda c, rl, i: (c,)), ('class,i', lambda c, rl, i: (c, i)),
              ('class,cov-left/right/both', lambda c, rl, i: (c,) + feat(rl, i)),
              ('class,i,cov', lambda c, rl, i: (c, i) + feat(rl, i)),
              ('class,uncovered vertex', lambda c, rl, i: (c, feat(rl, i) == ('cl0', 'cr0', 'both0') or (feat(rl,i)[0] in ('cl0','clNone') and feat(rl,i)[1] in ('cr0','crNone'))))]:
    # recompute with row keys
    pass
def det(name, f, tgt):
    g = collections.defaultdict(collections.Counter)
    for rs in rows:
        for i in range(1, N + 1):
            g[f(rs, i)][tgt(rs, i)] += 1
    tot = sum(sum(h.values()) for h in g.values()); ok = sum(max(h.values()) for h in g.values())
    # leave-one-out: predict each cell by the modal image of its group minus itself; singleton group = miss (no prediction)
    loo = 0; loo_fb = 0
    cg = collections.defaultdict(collections.Counter)
    for rs in rows:
        for i in range(1, N + 1): cg[C[rs]][tgt(rs, i)] += 1
    for rs in rows:
        for i in range(1, N + 1):
            h = g[f(rs, i)]; t = tgt(rs, i)
            hh = collections.Counter(h); hh[t] -= 1
            hh = +hh
            if hh:
                m = max(hh.items(), key=lambda kv: (kv[1], kv[0]))[0]
                loo += (m == t); loo_fb += (m == t)
            else:
                ch = collections.Counter(cg[C[rs]]); ch[t] -= 1; ch = +ch
                loo_fb += (max(ch.items(), key=lambda kv: (kv[1], kv[0]))[0] == t) if ch else 0
    print('   leave-one-out %-9s: accuracy singleton=miss %.4f ; singleton->class-modal fallback %.4f ; singleton cells=%d' % (tgt.__name__, loo / tot, loo_fb / tot, sum(1 for h in g.values() if sum(h.values()) == 1)))
    print('%-40s -> %-9s groups=%4d modal share=%.4f pure=%d' % (name, tgt.__name__ if hasattr(tgt,'__name__') else '', len(g), ok / tot, sum(1 for h in g.values() if len(h) == 1)))
def tcl(rs, i): return img[rs, i][0]
def tky(rs, i): return img[rs, i][1]
unc = lambda rs, i: feat(R[rs], i)
fs = [('class', lambda rs, i: C[rs]), ('class,i', lambda rs, i: (C[rs], i)), ('class,cov', lambda rs, i: (C[rs], unc(rs, i))),
      ('class,i,cov', lambda rs, i: (C[rs], i, unc(rs, i)))]
for nm, f in fs:
    det(nm, f, tcl); det(nm, f, tky)
# free-stretch filter: vertex i at distance >= d from every nonzero relation's support [v, v+r]
def dist(rl, i):
    ds = [0 if v <= i <= v + r else min(abs(i - v), abs(i - v - r)) for v, r in enumerate(rl, 1) if r]
    return min(ds) if ds else 99
print('\nTABLE 4: free-stretch filter. Cells (row,i) with dist(i, every relation support) >= d; per class: #image classes among filtered cells')
for d in (1, 2, 3):
    cells = [(rs, i) for rs in rows for i in range(1, N + 1) if dist(R[rs], i) >= d]
    g = collections.defaultdict(set)
    for rs, i in cells: g[C[rs]].add(img[rs, i][0])
    gk = collections.defaultdict(set)
    for rs, i in cells: gk[C[rs]].add(img[rs, i][1])
    print('d>=%d: cells=%d classes with cells=%d/%d; classes with one image class=%d, one image key=%d; mean #images/class %.2f' % (
        d, len(cells), len(g), len(cls), sum(1 for v in g.values() if len(v) == 1), sum(1 for v in gk.values() if len(v) == 1), sum(len(v) for v in g.values()) / max(1, len(g))))
    for c in cls:
        if c in g: print('    %-22s cells=%4d images=%s' % (c, sum(1 for rs, i in cells if C[rs] == c), sorted(L[x] for x in g[c])))
# pair compatibility under filter: same-class pairs of rows with filtered i, i' ; fraction with equal image class (any filtered i,i')
print('\nTABLE 5: same-class pairs, both rows having a filtered vertex (d>=2): all (i,i\') choices agree / some agree')
for d in (1, 2, 3):
    allag = some = pairs = 0
    for c, rs_ in by.items():
        S = {}
        for rs in rs_:
            s = {img[rs, i][0] for i in range(1, N + 1) if dist(R[rs], i) >= d}
            if s: S[rs] = s
        L_ = list(S.values()); 
        cnt = collections.Counter(frozenset(s) for s in L_)
        keys = list(cnt)
        for a in range(len(keys)):
            for b in range(a, len(keys)):
                m = cnt[keys[a]] * (cnt[keys[a]] - 1) // 2 if a == b else cnt[keys[a]] * cnt[keys[b]]
                if not m: continue
                pairs += m
                if keys[a] & keys[b]: some += m
                if len(keys[a]) == 1 and keys[a] == keys[b]: allag += m
    print('d>=%d: pairs=%d  all-choices-agree (singletons equal)=%d (%.3f)  some-choice-agrees=%d (%.3f)' % (d, pairs, allag, allag / pairs, some, some / pairs))
