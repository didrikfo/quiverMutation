"""T4 round 098 response: (1) r_v at step 1 for the 10 width-5 A rules (80 matches) and the 14 B rules, with the rule's window pattern;
(2) the 6 exceptions (B by replay, first step LNA): rule, matches, r_v at step 1, failing step and reason, lengths w+1..w+3;
(3) count of B-by-replay rules whose first step leaves, at w+1..w+3 (E-262's 304)."""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
def rv_of(L, r, s0):
    rr = r if s0 > 0 else rev(L, r); v = abs(s0) if s0 > 0 else L + 1 - abs(s0)
    return v, (rr[v - 1] if 1 <= v <= L - 2 else 'end')
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None]
print('== 1. width 5: r_v at step 1 per rule (lengths 6..8)')
for d in [x for x in flo if x[0] == 5]:
    c = collections.Counter(); cls = collections.Counter()
    for L in range(6, 9):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                o = replay(L, r, d, ws)
                if o == 'noapply': continue
                s0 = lm.applyAt(L, r, d, ws)[1][0]
                c[(rv_of(L, r, s0)[1], step(L, r, s0)[0])] += 1; cls['A' if o == 'agree' else 'B'] += 1
    print(d[1:3], 'seq', d[3], dict(cls), 'r_v,status', dict(c))
print('== 2. exceptions, lengths w+1..w+3')
exc = []; b_first = collections.Counter(); nb = 0
for d in flo:
    full = collections.Counter(); first = collections.Counter(); rows = []
    for L in range(d[0] + 1, d[0] + 4):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                o = replay(L, r, d, ws)
                if o == 'noapply': continue
                s0 = lm.applyAt(L, r, d, ws)[1][0]
                st, _ = step(L, r, s0); first[st] += 1; full['A' if o == 'agree' else 'B'] += 1
                rows.append((L, r, ws, s0, o, st))
    if full['B'] and not full['A']:
        nb += 1
        if first['notLNA']: b_first['leaves'] += 1
        else:
            b_first['first step LNA'] += 1; exc.append((d, rows))
print('B-by-replay rules at w+1..w+3:', nb, dict(b_first))
for d, rows in exc:
    print('RULE', d[:3], 'seq', d[3], 'matches', len(rows))
    for L, r, ws, s0, o, st in rows[:3]:
        v, rv = rv_of(L, r, s0)
        # state after step 1, and reason at step 2
        ap = lm.applyAt(L, r, d, ws)[1]; _, r1 = step(L, r, s0)
        order = list(range(1, L + 1)); pos = order.index(abs(s0)) + 1
        a = pos - 1 if s0 > 0 else pos - 2; order[a], order[a + 1] = order[a + 1], order[a]
        pos2 = order.index(abs(ap[1])) + 1; s2 = pos2 if ap[1] > 0 else -pos2
        n_ = L; rr1 = r1 if s2 > 0 else rev(L, r1); v2 = pos2 if s2 > 0 else L + 1 - pos2
        r2 = rr1[v2 - 1] if 1 <= v2 <= L - 2 else 'end'
        print('  L', L, 'row', r, 'seq', ap, 'step1 s', s0, 'r_v', rv, '(', st, ') row after', r1, 'step2 s', s2, 'r_v2', r2, 'outcome', o)
