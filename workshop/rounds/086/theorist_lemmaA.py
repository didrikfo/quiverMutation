"""Round 086 theorist, T5: machine check of Lemma A (monomial case) on quivers with n <= 5 vertices.
Statement checked: T = path algebra of a graded quiver Q (n vertices, any arrows incl. loops/parallel, degrees 0..2),
X1, Y disjoint arrow subsets, THETA = one theta_{xy} (deg |x|+|y|+1, weight 2) for each composable pair x in X1, y in Y
(x then y).  D theta = xy, D = 0 on arrows, D a derivation (Koszul sign (-1)^{deg of prefix}).
Then H(T<theta>) in each (weight, s, t, degree) = number of words in Q avoiding every bad pair xy (= basis of T/(xy)).
Also checks the tiling fact used in the direct proof: bad sites of a word are disjoint, and D keeps the tile skeleton.
Usage: .venv/bin/python workshop/rounds/086/theorist_lemmaA.py [trials_per_n] [WMAX] [seed]
"""
import sys, random
P = 1000003
def rank_mod(rows):
    piv = {}; r = 0
    for row in rows:
        row = {c: v % P for c, v in row.items() if v % P}
        while row:
            c = min(row)
            if c in piv:
                f = row[c]
                for cc, vv in piv[c].items():
                    nv = (row.get(cc, 0) - f * vv) % P
                    if nv: row[cc] = nv
                    else: row.pop(cc, None)
            else:
                inv = pow(row[c], P - 2, P)
                piv[c] = {cc: vv * inv % P for cc, vv in row.items()}; r += 1; break
    return r

def check(arrows, X1, Y, WMAX):
    """arrows: list of (s,t,deg). letters: ('a',k) arrows, ('t',x,y) thetas. returns (cells, bad)"""
    pairs = [(x, y) for x in X1 for y in Y if arrows[x][1] == arrows[y][0]]
    pset = set(pairs)
    L = [('a', k) for k in range(len(arrows))] + [('t', x, y) for x, y in pairs]
    def src(l): return arrows[l[1]][0] if l[0] == 'a' else arrows[l[1]][0]
    def tgt(l): return arrows[l[1]][1] if l[0] == 'a' else arrows[l[2]][1]
    def dg(l): return arrows[l[1]][2] if l[0] == 'a' else arrows[l[1]][2] + arrows[l[2]][2] + 1
    def wt(l): return 1 if l[0] == 'a' else 2
    words = {}   # weight -> list of words
    words[0] = [()]
    byw = {0: [((), None)]}
    allw = {0: [()]}
    for w in range(1, WMAX + 1):
        out = []
        for l in L:
            if wt(l) <= w:
                for u in allw[w - wt(l)]:
                    # append l at the end
                    if u == () or tgt(u[-1]) == src(l): out.append(u + (l,))
        allw[w] = out
    cells = 0; bad = 0
    for w in range(1, WMAX + 1):
        grp = {}
        for u in allw[w]:
            if u == (): continue
            grp.setdefault((src(u[0]), tgt(u[-1])), []).append(u)
        for key, ws in grp.items():
            bydeg = {}
            for u in ws: bydeg.setdefault(sum(dg(l) for l in u), []).append(u)
            idx = {d: {u: n for n, u in enumerate(l)} for d, l in bydeg.items()}
            rk = {}
            for d, l in bydeg.items():
                rows = []
                for u in l:
                    b = {}; pre = 0
                    for k, l_ in enumerate(u):
                        if l_[0] == 't':
                            nu = u[:k] + (('a', l_[1]), ('a', l_[2])) + u[k + 1:]
                            b[nu] = b.get(nu, 0) + (-1) ** pre
                        pre += dg(l_)
                    tgtidx = idx.get(d - 1)
                    row = {}
                    for nu, c in b.items():
                        assert tgtidx is not None and nu in tgtidx, 'D leaves cell'
                        row[tgtidx[nu]] = c
                    rows.append(row)
                rk[d] = rank_mod(rows)
            # normal words: all letters arrows, no bad pair
            for d, l in bydeg.items():
                n = len(l)
                h = n - rk.get(d, 0) - rk.get(d + 1, 0)
                normal = sum(1 for u in l if all(x[0] == 'a' for x in u) and
                             not any((u[k][1], u[k + 1][1]) in pset for k in range(len(u) - 1)))
                cells += 1
                if h != normal: bad += 1
    # tiling check on all words: bad sites are disjoint, D preserves site positions
    tilebad = 0
    for w in allw:
        for u in allw[w]:
            sites = []
            for k, l in enumerate(u):
                if l[0] == 't': sites.append((k, k + 1))
            # in the expanded (all-arrow) form positions; recompute on expanded word
            e = []
            for l in u:
                e += [l] if l[0] == 'a' else [('a', l[1]), ('a', l[2])]
            ps = [k for k in range(len(e) - 1) if (e[k][1], e[k + 1][1]) in pset]
            if any(ps[i + 1] - ps[i] < 2 for i in range(len(ps) - 1)): tilebad += 1
    return cells, bad, tilebad

def rand_instance(rng, n):
    na = rng.randint(2, 6)
    arrows = [(rng.randrange(n), rng.randrange(n), rng.choice([0, 0, 1, 2])) for _ in range(na)]
    idx = list(range(na)); rng.shuffle(idx)
    c = rng.randint(1, na - 1) if na > 1 else 1
    X1 = set(idx[:c]); Y = set(idx[c:]) if rng.random() < .7 else set(idx[c:c + 1])
    # also allow arrows in neither
    if rng.random() < .5 and len(Y) > 1: Y.pop()
    return arrows, sorted(X1), sorted(Y)

if __name__ == '__main__':
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    tot = [0, 0, 0, 0]
    for n in range(1, 6):
        c = b = tb = npairs = 0
        for _ in range(T):
            arrows, X1, Y = rand_instance(rng, n)
            assert not set(X1) & set(Y)
            cc, bb, t = check(arrows, X1, Y, W)
            c += cc; b += bb; tb += t
            npairs += sum(1 for x in X1 for y in Y if arrows[x][1] == arrows[y][0])
        print(f'n={n}: instances={T} (x,y) pairs={npairs} cells={c} mismatches={b} overlapping-sites={tb}')
        tot[0] += T; tot[1] += c; tot[2] += b; tot[3] += tb
    print('TOTAL instances, cells, mismatches, overlaps:', tot)
    # negative control: X1 and Y overlap (an arrow in both) -> Lemma A should fail (mismatches / overlapping sites)
    cb = ct = ci = 0
    for _ in range(4 * T):
        arrows, X1, Y = rand_instance(rng, rng.randint(1, 3))
        z = rng.choice(X1); Y = sorted(set(Y) | {z})
        if not any(arrows[x][1] == arrows[y][0] for x in X1 for y in Y) or not (set(X1) & set(Y)): continue
        try:
            cc, bb, t = check(arrows, X1, Y, W)
        except AssertionError:
            bb = t = 1
        ci += 1; cb += (bb > 0); ct += (t > 0)
    print(f'CONTROL (X1 and Y share an arrow): instances={ci}, with mismatch or D-leaves-cell={cb}, with overlapping sites={ct}')
