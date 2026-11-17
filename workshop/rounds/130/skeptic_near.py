"""Skeptic r130 response, item 3: the sharp locality audit.  Same matches as skeptic_h029.py (class-B floating rules,
lengths w+1..w+1+EXTRA), split by whether an OUTSIDE relation (path j -> j+a not inside the window ws..ws+w) has an endpoint
within 1 vertex of |s0|: min(|v-j|, |v-j-a|) <= 1, v = |s0|.  Each such match is compared with the rule's isolated (L = w+1)
value.  Also tallies, per rule, the nearest outside-relation endpoint to |s0| ever seen.
usage: skeptic_near.py W1,W2,... [EXTRA=3]"""
import sys, collections, time
MYARGV = list(sys.argv)
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import lm, nakayama, mutation, rowsOf, step, rels

W = [int(a) for a in MYARGV[1].split(',')]
EXTRA = int(MYARGV[2]) if len(MYARGV) > 2 else 3

for w in W:
    t0 = time.time()
    flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == w]
    data = collections.defaultdict(list)          # rule -> [(L, dist, delta)]
    for i, d in enumerate(flo):
        for L in range(w + 1, w + 2 + EXTRA):
            for r in rowsOf(L):
                for ws in lm.windowStartsFor(L, d):
                    if not lm.matchesAt(L, r, d, ws): continue
                    ap = lm.applyAt(L, r, d, ws)
                    if ap is None: continue
                    s0 = ap[1][0]
                    st, _ = step(L, r, s0)
                    if st != 'notLNA': continue
                    v = abs(s0)
                    out = [(j, a) for j, a in rels(r) if j < ws or j + a > ws + w]
                    dist = min((min(abs(v - j), abs(v - j - a)) for j, a in out), default=None)
                    alg = nakayama.LinearNakayamaAlgebra(L, r)
                    m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                    y = (len(list(m.arrows())) - len(list(alg.arrows())), len(m.rels) - len(alg.rels))
                    data[i].append((L, dist, y))
    near = [];  nearrules = 0; closest = collections.Counter()
    for v in data.values():
        iso = set(x[2] for x in v if x[0] == w + 1)
        nv = [(x, x[2] in iso) for x in v if x[1] is not None and x[1] <= 1]
        near += nv; nearrules += bool(nv)
        ds = [x[1] for x in v if x[1] is not None]
        closest[min(ds) if ds else 'none'] += 1
    print('width %d lengths %d..%d: class-B rules %d, matches %d' % (w, w + 1, w + 1 + EXTRA, len(data), sum(map(len, data.values()))))
    print('  outside relation within 1 of |s0|: %d matches on %d rules, values %s, differ from isolated %d'
          % (len(near), nearrules, dict(collections.Counter(x[2] for x, _ in near)), sum(not ok for _, ok in near)))
    print('  rules by closest outside-relation endpoint to |s0|:', dict(sorted(closest.items(), key=str)))
    print('  time %.0f s' % (time.time() - t0)); sys.stdout.flush()
