"""T4 round 095: does E-187's one-step row rule (and its E-205 mirror) alone reproduce every floating H-020 rule?
Row model: LNA on n vertices = row r[1..n-2] (r_j = length of relation starting at arrow j, 0 = none).
Right mutation at v (E-187): gate iff v<=n-1 and no relation ends at vertex v+1 (j+r_j=v+1);
 result is an LNA iff v=n-1 or r_v=2, then r_v:=0 (v<=n-2), r_{v-1}:=2 (v>=2).
Left mutation at v = right mutation at n+1-v of the reversed row (E-205 via _opposite); reversal of the
 relation (start j, len a) -> (start n-j-a+1, len a).
The engine keeps labels: a mutation swaps two adjacent line positions (right: v,v+1; left: v-1,v) and later steps name labels.
Usage: theorist_rowmodel.py check   (row model vs repo engine, n<=8, all rows, all v)
       theorist_rowmodel.py rules LMAX  (every floating rule: replay in row model vs applyAt prediction)"""
import sys, collections
from quivermutation import lnaMoves as lm, nakayama, mutation

def rels(r): return [(j + 1, a) for j, a in enumerate(r) if a]
def rev(n, r):
    out = [0] * (n - 2)
    for j, a in rels(r): out[n - j - a + 1 - 1] = a
    return out
def right(n, r, v):
    """returns ('gate',None) if gate fails, ('notLNA',None), or ('ok',newrow)"""
    if v > n - 1 or v < 1: return 'gate', None
    if any(j + a == v + 1 for j, a in rels(r)): return 'gate', None
    if not (v == n - 1 or r[v - 1] == 2): return 'notLNA', None
    new = list(r)
    if v <= n - 2: new[v - 1] = 0
    if v >= 2: new[v - 2] = 2
    return 'ok', new
def step(n, r, s):
    if s > 0: return right(n, r, s)
    st, new = right(n, rev(n, r), n + 1 + s)   # left at v=-s  <-> right at n+1-v on reversed
    return (st, rev(n, new) if new else None)

def allrows(n):
    out = []
    def rec(j, cur):
        if j == n - 2:
            if lm.isAdmissible(n, cur): out.append(list(cur))
            return
        for a in [0] + list(range(2, n - j)):
            cur.append(a)
            rs = rels(cur)
            if all(e[0] + e[1] < l[0] + l[1] for e, l in zip(rs, rs[1:])): rec(j + 1, cur)
            cur.pop()
    rec(0, []); return out

def check():
    bad = collections.Counter(); tot = 0
    for n in range(4, 9):
        for r in allrows(n):
            alg = nakayama.LinearNakayamaAlgebra(n, r)
            for s in list(range(1, n + 1)) + [-v for v in range(1, n + 1)]:
                target = alg if s > 0 else lm._quiet(lm.pathAlgebra.dualPathAlgebra, alg)
                gate = bool(lm._quiet(mutation.mutationIsPossibleAtVertex, target, abs(s)))
                st, new = step(n, r, s); tot += 1
                if gate != (st != 'gate'): bad['gate', s > 0] += 1; continue
                if not gate: continue
                m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s])
                act = lm.asRelLengths(m, n)
                if (act is not None) != (st == 'ok'): bad['lna', s > 0] += 1
                elif act is not None and list(act) != new: bad['image', s > 0] += 1
    print('row model vs engine: cases', tot, 'violations', dict(bad) or 'none')

_cache = {}
def rowsOf(n):
    if n not in _cache: _cache[n] = allrows(n)
    return _cache[n]

def replay(L, r, d, ws):
    """'agree' | 'DISAGREE' | 'noapply' | (stepIndex, 'gate'/'notLNA')."""
    ap = lm.applyAt(L, r, d, ws)
    if ap is None: return 'noapply'
    pred, seq = ap; cur = r; order = list(range(1, L + 1))
    for k, s in enumerate(seq):
        pos = order.index(abs(s)) + 1       # position of the engine label in the line
        st, nx = step(L, cur, pos if s > 0 else -pos)
        if st != 'ok': return (k, st)
        cur = nx
        a = pos - 1 if s > 0 else pos - 2   # E-187/E-205 relabelling: right swaps positions v,v+1; left v-1,v
        order[a], order[a + 1] = order[a + 1], order[a]
    return 'agree' if cur == pred else 'DISAGREE'

def rules(widths, extra):
    """Every floating rule of the given widths, lengths w+1..w+extra, all matches: replay in the row model."""
    flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] in widths]
    res = {}
    for d in flo:
        w = d[0]; c = collections.Counter()
        for L in range(w + 1, w + extra + 1):
            for r in rowsOf(L):
                for ws in lm.windowStartsFor(L, d):
                    if lm.matchesAt(L, r, d, ws):
                        out = replay(L, r, d, ws)
                        c[out if isinstance(out, str) else out[1] + '@step%d' % out[0]] += 1
        res[d] = c
    return res

if __name__ == '__main__':
    if sys.argv[1] == 'check': check()
    else:
        widths = [int(x) for x in sys.argv[2].split(',')]; extra = int(sys.argv[3])
        res = rules(widths, extra)
        cls = collections.Counter(); tot = collections.Counter()
        for d, c in res.items():
            bad = sum(v for k, v in c.items() if k != 'agree')
            kind = 'A all-LNA' if c['agree'] and not bad else ('B never-LNA' if not c['agree'] else 'C mixed')
            if not c: kind = 'D no match'
            cls[(d[0], kind)] += 1; tot.update(c)
            if c['DISAGREE']: print('DISAGREE', d, dict(c))
        print('widths', widths, 'lengths w+1..w+%d' % extra)
        for k in sorted(cls): print(k, cls[k])
        print('match totals', dict(tot))
