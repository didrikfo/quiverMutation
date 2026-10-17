"""T7 (round 070): left mutation of an LNA row, checked directly (not via the dual rule).
Usage: theorist_leftcheck.py NMIN NMAX.  Uses mutation.quiverMutationAtVertices(alg,[-v]) (repo left mutation,
= procedure.mutateLeftAtVertex) and the gate on the dual algebra, for every row r and v = 1..n.
Predicted (mirror of E-187): relation (j, j+r_j) has r_j arrows.
 G: gate <=> v >= 2 and no relation STARTS at v-1.
 L: (gate) LNA <=> v = 2 or r_{v-2} = 2 (relation (v-2, v)); image: r_{v-2} := 0 (v>=3), r_{v-1} := 2 (v <= n-1).
 C: per-vertex LNA count = C_{v-2} C_{n-v}.   F: otherwise non-LNA (arrow count n-1 + out/in-degree 2 or n arrows)."""
import sys, collections
from math import comb
sys.path.insert(0, 'workshop/rounds/062')
from quivermutation import lnaMoves as lm, mutation, nakayama, pathAlgebra
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
cat = lambda k: comb(2 * k, k) // (k + 1)
bad = collections.Counter(); cnt = collections.Counter(); pv = collections.Counter()
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for r in allLNAs(n):
        starts = {j + 1 for j in range(n - 2) if r[j]}
        alg = nakayama.LinearNakayamaAlgebra(n, r)
        dual = lm._quiet(pathAlgebra.dualPathAlgebra, alg)
        for v in range(1, n + 1):
            gate = bool(lm._quiet(mutation.mutationIsPossibleAtVertex, dual, v))
            if gate != (v >= 2 and (v - 1) not in starts): bad['gate', n] += 1
            if not gate: continue
            nx_ = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [-v])
            new = lm.asRelLengths(nx_, n)
            want = (v == 2 or (v >= 3 and r[v - 3] == 2))
            cnt[n, 'LNA' if new is not None else 'not'] += 1
            if (new is not None) != want: bad['lna', n] += 1
            if new is not None:
                pv[n, v] += 1
                img = list(r)
                if v >= 3: img[v - 3] = 0
                if v <= n - 1: img[v - 2] = 2
                if img != new: bad['image', n] += 1
            else:
                q = nx_.quiver; E = len(q.edges)
                od = max(max(d for _, d in q.out_degree()), max(d for _, d in q.in_degree()))
                if not ((E == n - 1 or E == n) and od == 2): bad['F', n] += 1
                if want: bad['unexpected-fail', n] += 1
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for v in range(2, n + 1):
        if pv[n, v] != cat(v - 2) * cat(n - v): bad['C', n] += 1
    print('C per-vertex n =', n, [pv[n, v] for v in range(2, n + 1)], 'total', sum(pv[n, v] for v in range(2, n + 1)), 'C_{n-1} =', cat(n - 1))
for k in sorted(cnt): print(k, cnt[k])
print("violations:", dict(bad) or "none")
