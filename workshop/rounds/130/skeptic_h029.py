"""Skeptic r130, T4 / H-029 at widths 8-9: is the first-step (arrows, relations) delta of a class-B floating rule a function
of the rule alone?  Pre-registered (before running widths 8-9): YES = every notLNA rule has one delta over all its matches at
lengths w+1..w+3 and all window positions; and the delta is predicted by the rule's ISOLATED value (length w+1, window = whole
line, no context) wherever that exists.  One rule with two deltas, or one wrong isolated prediction, refutes H-029 at that width.
Locality audit: matchesAt forbids relations straddling the window, so the only context that can reach the first step is a
relation that touches a window END VERTEX (ends at vertex ws or starts at vertex ws+w).  Matches are split by that.
Controls (power): (a) labels shuffled over all matches of a width, rule match-counts kept -> how many rules stay constant;
(b) one planted wrong delta per rule in a random half of the rules -> how many caught.
usage: skeptic_h029.py W1,W2,... [EXTRA=2] [NSHUF=200]     (lengths w+1 .. w+1+EXTRA)"""
import sys, collections, random, time
MYARGV = list(sys.argv)
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import lm, nakayama, mutation, rowsOf, step, rels

W = [int(a) for a in MYARGV[1].split(',')]
EXTRA = int(MYARGV[2]) if len(MYARGV) > 2 else 2
NSHUF = int(MYARGV[3]) if len(MYARGV) > 3 else 200
rng = random.Random(130)


def touches(L, r, ws, w):
    """an outside relation shares an end vertex with the window (window vertices ws..ws+w)"""
    for j, a in rels(r):
        if j + a == ws or j == ws + w: return True
    return False


for w in W:
    t0 = time.time()
    flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == w]
    data = collections.defaultdict(list)          # rule index -> [(L, ws, touch, delta)]
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
                    alg = nakayama.LinearNakayamaAlgebra(L, r)
                    m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                    y = (len(list(m.arrows())) - len(list(alg.arrows())), len(m.rels) - len(alg.rels))
                    data[i].append((L, ws, touches(L, r, ws, w), y))
    nrule = len(data); nmatch = sum(len(v) for v in data.values())
    const = sum(len(set(x[3] for x in v)) == 1 for v in data.values())
    multi = [(flo[i][1:4], collections.Counter((x[0], x[2], x[3]) for x in v)) for i, v in data.items() if len(set(x[3] for x in v)) > 1]
    # isolated prediction: value at L = w+1 (no context at all) predicts every longer match
    iso_n = iso_ok = iso_rules = 0
    for v in data.values():
        base = set(x[3] for x in v if x[0] == w + 1)
        if len(base) != 1: continue
        iso_rules += 1
        for x in v:
            if x[0] > w + 1: iso_n += 1; iso_ok += x[3] in base
    # length transfer as in E-310: deltas seen at lengths < L predict L
    tr_n = tr_ok = 0
    for v in data.values():
        Ls = sorted(set(x[0] for x in v))
        for L in Ls[1:]:
            hist = set(x[3] for x in v if x[0] < L)
            if len(hist) == 1: tr_n += 1; tr_ok += set(x[3] for x in v if x[0] == L) <= hist
    tch = [x for v in data.values() for x in v if x[2]]
    tch_rules = sum(any(x[2] for x in v) and any(not x[2] for x in v) for v in data.values())
    mrule = collections.Counter(len(v) for v in data.values())
    print('width %d lengths %d..%d: floating %d, notLNA-first rules %d, matches %d (per rule %s); one delta %d, several %d'
          % (w, w + 1, w + 1 + EXTRA, len(flo), nrule, nmatch, dict(sorted(mrule.items())), const, len(multi)))
    print('  delta values', dict(collections.Counter(x[3] for v in data.values() for x in v)))
    print('  isolated (L=w+1) value predicts longer matches: %d / %d (rules with an isolated value %d)' % (iso_ok, iso_n, iso_rules))
    print('  length transfer: %d / %d' % (tr_ok, tr_n))
    print('  matches with context on a window end vertex: %d; rules seen both with and without it: %d' % (len(tch), tch_rules))
    for rule, c in multi[:8]: print('  MULTI', rule, dict(c))
    # controls
    labels = [x[3] for v in data.values() for x in v]; sizes = [len(v) for v in data.values()]
    nonsingle = sum(s >= 2 for s in sizes)
    shuf_const = []
    for _ in range(NSHUF):
        rng.shuffle(labels); k = 0; c = 0
        for s in sizes:
            c += len(set(labels[k:k + s])) == 1; k += s
        shuf_const.append(c)
    print('  control (a) shuffled labels, %d shuffles: constant rules max %d, mean %.1f (observed %d; rules with >= 2 matches %d)'
          % (NSHUF, max(shuf_const) if shuf_const else -1, sum(shuf_const) / max(1, len(shuf_const)), const, nonsingle))
    vals = sorted(set(labels)) or [(1, 1)]
    planted = caught = 0
    for v in data.values():
        if len(v) < 2 or rng.random() < 0.5: continue
        ys = [x[3] for x in v]; j = rng.randrange(len(ys))
        ys[j] = (ys[j][0], ys[j][1] + 1)        # planted wrong delta in one match
        planted += 1; caught += len(set(ys)) > 1
    print('  control (b) planted one wrong delta: caught %d / %d' % (caught, planted))
    print('  time %.0f s' % (time.time() - t0)); sys.stdout.flush()
