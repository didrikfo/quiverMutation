"""Skeptic r106 (T4/H-020): is there a one-step rule for the middle (non-LNA) algebra of the 14 width-5 class-B rules that a control catches?
Candidate rule R (pre-registered): after the first mutation at v with r_v = k >= 3 (gate ok), the mutated algebra has
  arrows = (n-1) + (k-2)  and  the row-model count of surviving relations is unchanged;  i.e. arrows depend on (n,k) only.
Test: engine arrow count of the middle algebra over all matches of the 14 B rules (lengths 6..8) vs R.
Controls: (a) planted wrong rule R' (arrows = n-1+1 for all k) must fail where k=4 occurs; (b) the 10 A rules (r_v=2, LNA) must give arrows = n-1 (rule with k=2)."""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == 5]

def feat(m): return (len(list(m.arrows())), len(list(m.rels)))
res = collections.defaultdict(collections.Counter); cls = {}
for d in flo:
    for L in range(6, 9):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                s0 = lm.applyAt(L, r, d, ws)[1][0]
                st, _ = step(L, r, s0)
                if st == 'gate': continue
                rr = r if s0 > 0 else rev(L, r); v = abs(s0) if s0 > 0 else L + 1 - abs(s0)
                k = rr[v-1] if 1 <= v <= L-2 else 'end'
                alg = nakayama.LinearNakayamaAlgebra(L, r)
                m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                a0 = len(list(alg.arrows())); r0 = len(list(alg.rels))
                f = feat(m)
                res[(d[1:3], k)][(f[0]-a0, f[1]-r0)] += 1
for key, c in sorted(res.items(), key=lambda x: str(x)): print(key, dict(c))

# ---- stage 2: does a local feature tuple DETERMINE the (d_arrows, d_rels) change? fit on L=6,7, predict L=8 (held-out length).
# Pre-registered feature F = (k, r_{v-1}>0, r_{v+1}>0) of the row in mutation direction.  Controls: F0 = (k,) only; Fbad = (k, L%2, hash of row) (planted junk).
rows_data = []
for d in flo:
    for L in range(6, 9):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                s0 = lm.applyAt(L, r, d, ws)[1][0]
                st, _ = step(L, r, s0)
                if st != 'notLNA': continue
                rr = r if s0 > 0 else rev(L, r); v = abs(s0) if s0 > 0 else L + 1 - abs(s0)
                g = lambda j: rr[j-1] if 1 <= j <= L-2 else 0
                alg = nakayama.LinearNakayamaAlgebra(L, r)
                m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                y = (len(list(m.arrows())) - len(list(alg.arrows())), len(m.rels) - len(alg.rels))
                rows_data.append((L, d, tuple(rr), v, g(v), g(v-1), g(v+1), y))
print('notLNA matches', len(rows_data))
def test(name, F):
    tab = collections.defaultdict(collections.Counter)
    for x in rows_data:
        if x[0] <= 7: tab[F(x)][x[-1]] += 1
    hit = miss = unseen = 0
    for x in rows_data:
        if x[0] == 8:
            c = tab.get(F(x))
            if not c: unseen += 1
            elif c.most_common(1)[0][0] == x[-1]: hit += 1
            else: miss += 1
    amb = sum(1 for c in tab.values() if len(c) > 1)
    print(f'{name}: train cells {len(tab)} ambiguous {amb}; L=8 hit {hit} miss {miss} unseen {unseen}')
test('F  (k, r_{v-1}>0, r_{v+1}>0)', lambda x: (x[4], x[5] > 0, x[6] > 0))
test('F0 (k)', lambda x: (x[4],))
test('Fbad (k, row hash)  [planted junk]', lambda x: (x[4], hash(x[2])))
test('Fneg (k, r_{v-1}>0)', lambda x: (x[4], x[5] > 0))
test('Fpos (k, r_{v+1}>0)', lambda x: (x[4], x[6] > 0))
