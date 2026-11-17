"""Skeptic r130, T4 / H-029: H-029's falsifier reads 'build the rule table at widths <= 7 and predict the delta at widths 8-9'.
A width-8 rule is a different rule from every width-7 rule unless it is a PADDING of one (one free arrow added at the left or
right end of the window, pattern, image and offsets shifted).  This script finds, for every floating rule of width w, whether
it is a padding of a floating rule of width w-1 (recursively, down to the narrowest), and compares the first-step deltas
(computed as in skeptic_h029.py at lengths w+1..w+1+EXTRA).  Prints how many width-8/9 rules the <= 7 table can predict at all,
and how many it predicts right.  Control: the same prediction with deltas shuffled over the narrower rules.
usage: skeptic_h029pad.py WMAX [EXTRA=1]"""
import sys, collections, random
MYARGV = list(sys.argv)
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import lm, nakayama, mutation, rowsOf, step

WMAX = int(MYARGV[1]); EXTRA = int(MYARGV[2]) if len(MYARGV) > 2 else 1
flo = {w: [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == w] for w in range(3, WMAX + 1)}


def key(d): return (d[0], tuple(sorted(d[1])), tuple(sorted(d[2])), tuple(d[3]))
def padR(d): return (d[0] + 1, d[1], d[2], d[3])
def padL(d):
    sh = lambda rs: tuple((s + 1, a) for s, a in rs)
    return (d[0] + 1, sh(d[1]), sh(d[2]), tuple(v + 1 if v > 0 else v - 1 for v in d[3]))


def deltas(d):
    w = d[0]; out = set()
    for L in range(w + 1, w + 2 + EXTRA):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                ap = lm.applyAt(L, r, d, ws)
                if ap is None: continue
                s0 = ap[1][0]
                if step(L, r, s0)[0] != 'notLNA': return None if not out else out   # class A at this match
                alg = nakayama.LinearNakayamaAlgebra(L, r)
                m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                out.add((len(list(m.arrows())) - len(list(alg.arrows())), len(m.rels) - len(alg.rels)))
    return out


D = {}
for w in flo:
    for d in flo[w]:
        dd = deltas(d)
        if dd: D[key(d)] = dd
idx = {key(d): d for w in flo for d in flo[w]}
rng = random.Random(130)
for w in range(4, WMAX + 1):
    parents = {}
    for d in flo[w - 1]:
        for p in (padR(d), padL(d)): parents.setdefault(key(p), []).append(key(d))
    B = [key(d) for d in flo[w] if key(d) in D]
    padded = [k for k in B if k in parents]
    withB = [k for k in padded if any(p in D for p in parents[k])]
    ok = sum(any(D.get(p) == D[k] for p in parents[k]) for k in withB)
    allpad = sum(key(d) in parents for d in flo[w])
    # shuffled control: parent deltas permuted among the width-(w-1) class-B rules
    pk = [key(d) for d in flo[w - 1] if key(d) in D]; vals = [D[k] for k in pk]
    sh = []
    for _ in range(200):
        rng.shuffle(vals); P = dict(zip(pk, vals))
        sh.append(sum(any(P.get(p) == D[k] for p in parents[k] if p in P) for k in withB))
    print('width %d: floating %d, of them paddings of a width-%d floating rule %d; class-B %d, class-B paddings of a class-B rule %d,'
          ' delta equal to a parent %d; shuffled-parent control max %d mean %.1f'
          % (w, len(flo[w]), w - 1, allpad, len(B), len(withB), ok, max(sh) if sh else -1, sum(sh) / max(1, len(sh))))
    new = [idx[k][1:] for k in B if k not in parents]
    if new: print('   class-B rules that are not paddings (first 5 of %d):' % len(new), new[:5])


# Part 2: width families.  A rule's TEMPLATE replaces every relation length a >= 4 with a = w - c by ('W', c) and every offset
# v with |v| >= w - 1 by ('W', sign, w - |v|).  Rules of different widths with one template form a family; the <= w-1 table
# predicts a width-w rule when its template occurs at a smaller width with one delta.
def tmpl(d):
    w = d[0]
    L = lambda rs: tuple((s, ('W', w - a) if a >= 4 else a) for s, a in rs)
    o = tuple(('W', v > 0, w - abs(v)) if abs(v) >= w - 1 else v for v in d[3])
    return (L(d[1]), L(d[2]), o)


fam = collections.defaultdict(dict)
for k, dd in D.items(): fam[tmpl(idx[k])][k[0]] = dd
for w in range(6, WMAX + 1):
    B = [key(d) for d in flo[w] if key(d) in D]
    pred = [k for k in B if any(v for ww, v in fam[tmpl(idx[k])].items() if ww < w)]
    good = sum(all(v == D[k] for ww, v in fam[tmpl(idx[k])].items() if ww < w) for k in pred)
    # control: deltas of the smaller-width family members shuffled over all templates
    tl = [t for t in fam]; vals = [fam[t] for t in tl]
    sh = []
    for _ in range(200):
        perm = vals[:]; rng.shuffle(perm); F = dict(zip(tl, perm))
        sh.append(sum(any(ww < w for ww in F[tmpl(idx[k])]) and all(v == D[k] for ww, v in F[tmpl(idx[k])].items() if ww < w) for k in pred))
    print('width %d families: class-B %d, template seen at a smaller width %d, predicted right %d; shuffled-template control max %d mean %.1f'
          % (w, len(B), len(pred), good, max(sh), sum(sh) / len(sh)))
    bad = [(idx[k][1:], {ww: v for ww, v in fam[tmpl(idx[k])].items()}) for k in pred
           if not all(v == D[k] for ww, v in fam[tmpl(idx[k])].items() if ww < w)]
    for b in bad[:4]: print('   MISS', b)
