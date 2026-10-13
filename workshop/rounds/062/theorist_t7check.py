"""T7 (round 062): check the case analysis of the single-mutation LNA rule against the repo's procedure.
Usage: theorist_t7check.py NMIN NMAX.  For every LNA row r of length n and every v = 1..n, right mutation:
 G: gate  <=>  v <= n-1 and no relation ends at vertex v+1
 L: (gate) result is an LNA  <=>  v = n-1 or r_v = 2;  then row = r with (r_{v-1}, r_v) := (2, 0)  (v=1: r_1 := 0 only)
 F: (gate, not LNA) r_v = 0, v <= n-2: n-1 arrows and some vertex of out-degree 2;  r_v >= 3: n arrows (and out-degree 2 vertex)."""
import sys, collections
from quivermutation import lnaMoves as lm, mutation, nakayama
def allLNAs(n):
    out = []
    def rec(j, cur):
        if j == n - 2:
            if lm.isAdmissible(n, cur): out.append(list(cur))
            return
        for a in [0] + list(range(2, n - 1 - j + 1)):
            cur.append(a); rels = lm.relationsOf(cur)
            if all(e[0] + e[1] < l[0] + l[1] for e, l in zip(rels, rels[1:])): rec(j + 1, cur)
            cur.pop()
    rec(0, []); return out
bad = collections.Counter(); cnt = collections.Counter(); pv = collections.Counter()
from math import comb
cat = lambda k: comb(2 * k, k) // (k + 1)
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for r in allLNAs(n):
        ends = {j + 1 + r[j] for j in range(n - 2) if r[j]}
        alg = nakayama.LinearNakayamaAlgebra(n, r)
        for v in range(1, n + 1):
            gate = bool(lm._quiet(mutation.mutationIsPossibleAtVertex, alg, v))
            if gate != (v <= n - 1 and (v + 1) not in ends): bad['gate', n] += 1
            if not gate: continue
            nx_ = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [v])
            new = lm.asRelLengths(nx_, n)
            rv = r[v - 1] if v <= n - 2 else None
            want = (v == n - 1 or rv == 2)
            cnt[n, 'gated', 'LNA' if new is not None else 'not'] += 1
            if (new is not None) != want: bad['lna', n] += 1
            if new is not None: pv[n, v] += 1
            if new is not None:
                img = list(r)
                if v <= n - 2: img[v - 1] = 0
                if v >= 2: img[v - 2] = 2
                if img != new: bad['image', n] += 1
            else:
                q = nx_.quiver; E = len(q.edges); od = max(d for _, d in q.out_degree())
                if rv == 0 and not (E == n - 1 and od == 2): bad['F0', n] += 1
                if rv is not None and rv >= 3 and not (E == n and od == 2): bad['F3', n] += 1
                if rv is None or rv == 2: bad['unexpected-fail', n] += 1
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for v in range(1, n):  # (C) per vertex: LNA results = C_{v-1} C_{n-1-v}
        if pv[n, v] != cat(v - 1) * cat(n - 1 - v): bad['C', n] += 1
    print('C per-vertex n =', n, [pv[n, v] for v in range(1, n)])
for k in sorted(cnt): print(k, cnt[k])
print("violations:", dict(bad) or "none")
