"""T4 round 098: the 14 width-5 class-B floating rules vs the E-187 (right) / E-205 (left, mirror) one-step rule.
Criterion (fixed before running): for every match at lengths 6..8, (H) the row-model first step reports 'notLNA' iff the engine's
first mutation (same label) is not an LNA, and gate agrees; (R) the reason r_v (the row entry at the mutated position, 0 or >=3)
is constant per rule or splits only by end-distance; (S) per rule, which steps of the sequence are LNA in the engine."""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == 5]
Bs = []
for d in flo:
    c = collections.Counter()
    for L in range(6, 9):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if lm.matchesAt(L, r, d, ws):
                    o = replay(L, r, d, ws); c[o if isinstance(o, str) else o[1]] += 1
    if c and not c['agree']: Bs.append(d)
print('class B width 5:', len(Bs), 'of', len(flo))
tot = collections.Counter()
for d in Bs:
    seq = d[3]; rows = collections.Counter(); hit = collections.Counter(); n_m = 0
    for L in range(6, LMAX + 1):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                n_m += 1
                s0 = lm.applyAt(L, r, d, ws)[1][0]; pos = abs(s0)   # sequence shifted to the window start; first step: label = position
                st, new = step(L, r, s0)
                alg = nakayama.LinearNakayamaAlgebra(L, r)
                target = alg if s0 > 0 else lm._quiet(lm.pathAlgebra.dualPathAlgebra, alg)
                gate = bool(lm._quiet(mutation.mutationIsPossibleAtVertex, target, pos))
                m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0]) if gate else None
                eng = 'gate' if not gate else ('ok' if lm.asRelLengths(m, L) is not None else 'notLNA')
                hit[(st == eng)] += 1
                # reason: entry of the row at the mutated vertex (mirror for left)
                rr = r if s0 > 0 else rev(L, r); v = pos if s0 > 0 else L + 1 - pos
                ent = rr[v - 1] if 1 <= v <= L - 2 else 'end'
                dist = (v - 1 if v - 1 < 2 else 2, (L - 1 - v) if L - 1 - v < 2 else 2)
                rows[(st, ent if ent in (0, 'end') else ('>=3' if ent >= 3 else 2))] += 1
    tot.update(hit)
    print(d[1:], 'seq', seq, 'matches', n_m, 'E187==engine step1', dict(hit), 'status/entry', dict(rows))
print('TOTAL step-1 agreement', dict(tot))
