"""Round 097 skeptic: census of ALL (parent, vertex) pairs on random acyclic parents with the E-248 pattern (relation i ~> t with arrow i -> t),
split by what the search does: gate (procedure.isMutable), illegal-relation test, key guard, and (independent) tilting + End(T) dims.
Usage: .venv/bin/python workshop/rounds/097/skeptic_omega_census.py NMIN NMAX PARENTS [seed]"""
import sys, random, warnings
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/097')
warnings.simplefilter('ignore')
from collections import Counter
from skeptic_omega import *
from quivermutation import reduction, search as sr, arrowPaths as arrowPaths

def trunc_cartan(q, rels, L):
    """dims of e_a A e_b for a (possibly cyclic) bound quiver, paths of arrow length <= L, ideal = span of u r v with |u r v| <= L (truncated: an UPPER bound on the true dimension if the true algebra is finite dimensional with all paths of length > L in the ideal... not proved; compare L and L+2)."""
    V = sorted(q.nodes); arrows = ap.arrowsOf(q)
    out = {v: [a for a in arrows if a[0] == v] for v in V}
    paths = {(a, b): [] for a in V for b in V}
    cnt = 0
    def dfs(start, cur, path):
        nonlocal cnt
        paths[(start, cur)].append(tuple(path)); cnt += 1
        if cnt > 60000: raise OverflowError
        if len(path) == L: return
        for a in out[cur]: dfs(start, a[1], path + [a])
    for v in V: dfs(v, v, [])
    C = [[0] * len(V) for _ in V]
    for x, a in enumerate(V):
        for y, b in enumerate(V):
            P = paths[(a, b)]
            if not P: continue
            idx = {p: k for k, p in enumerate(P)}
            rows = []
            for r in rels:
                rs, rt = ap.source(r), ap.target(r)
                rl = max(len(p) for p in r)
                for u in paths[(a, rs)]:
                    for w in paths[(rt, b)]:
                        if len(u) + rl + len(w) > L: continue
                        rows.append({idx[u + p + w]: Fraction(c) for p, c in r.items()})
            C[x][y] = len(P) - rank(rows)
    return C

def deg2pat(q, rels, v):
    """Anick chain a b c from v with ab, bc zero relations (Ext^3(S_v, S_t) != 0) and an arrow v -> t."""
    zero = {next(iter(r)) for r in rels if len(r) == 1 and len(next(iter(r))) == 2}
    for a in ap.arrowsOutOf(q, v):
        for b in ap.arrowsOutOf(q, a[1]):
            if (a, b) not in zero: continue
            for c in ap.arrowsOutOf(q, b[1]):
                if (b, c) in zero and any(x[1] == c[1] for x in ap.arrowsOutOf(q, v)): return True
    return False

def census(nmin, nmax, per, seed, pat2=False):
    rng = random.Random(seed); cnt = Counter(); bad = []
    for n in range(nmin, nmax + 1):
        for _ in range(per):
            P = gen_parent(rng, n, pat2=pat2)
            if P is None: continue
            q = P.quiver; rels = procedure.relationsFrom(P)
            baseKey = sr._coxeterKeyOrNone(P)
            for v in sorted(q.nodes):
                if not ap.arrowsOutOf(q, v): continue
                pat = any(ap.source(r) == v and any(a[1] == ap.target(r) for a in ap.arrowsOutOf(q, v)) for r in rels)
                gate = procedure.isMutable(q, rels, v)
                pat3 = any(ap.source(r) == v and len(next(iter(r))) == 3 and any(a[1] == ap.target(r) for a in ap.arrowsOutOf(q, v)) for r in rels) or pat2
                lab = ['n=%d' % n, 'DEG2PAT' if deg2pat(q, rels, v) else 'no-deg2pat', 'PAT' if pat else 'nopat', 'gate' if gate else 'gate-refused']
                try:
                    nq, nr = procedure.mutateAtVertex(q, rels, v)
                except Exception as e:
                    cnt[tuple(lab + ['mutate raised ' + type(e).__name__])] += 1; continue
                ch = procedure.toPathAlgebra(nq, nr)
                illegal = any(arrowPaths.isIllegalRelation(ch.quiver, r) for r in procedure.relationsFrom(ch))
                ch = reduction.reducePathAlgebra(ch)
                key = sr._coxeterKeyOrNone(ch)
                kept = key is None or key == baseKey
                lab += ['illegal' if illegal else 'legalrel', 'keykept' if kept else 'keydrop', 'cyclic child' if not acyc(ch.quiver) else 'acyclic child']
                taken = gate and not illegal and kept
                lab.append('TAKEN' if taken else 'not taken')
                if not acyc(ch.quiver) and taken:
                    Vs, res = stepTest(P, v, -1); N = len(Vs); Hm = res[0]
                    try:
                        Cs = [trunc_cartan(ch.quiver, procedure.relationsFrom(ch), L) for L in (10, 12)]
                        ok = [all(c[a][b] == Hm[a][b] for a in range(N) for b in range(N)) or all(c[a][b] == Hm[b][a] for a in range(N) for b in range(N)) for c in Cs]
                        lab += ['trunc dims = End(T): L10 %s L12 %s' % tuple(ok), 'stable' if Cs[0] == Cs[1] else 'unstable']
                        if not all(ok): bad.append((n, v, rels, Cs, Hm))
                    except OverflowError: lab += ['path overflow']
                if acyc(ch.quiver):
                    Vs, res = stepTest(P, v, -1); N = len(Vs)
                    C = ap.cartanMatrix(ch.quiver, procedure.relationsFrom(ch), exact=True)
                    Hm = res[0]
                    eq = all(C[a][b] == Hm[a][b] for a in range(N) for b in range(N)) or all(C[a][b] == Hm[b][a] for a in range(N) for b in range(N))
                    J = sum(map(sum, res[-1])) + sum(map(sum, res[1]))
                    lab += ['End=child' if eq else 'MISMATCH', 'J=0' if J == 0 else 'J!=0']
                    if not eq and taken: bad.append((n, v, rels))
                cnt[tuple(lab)] += 1
    for k in sorted(cnt): print(cnt[k], k)
    print('taken mismatches:', len(bad))
if __name__ == '__main__':
    a = sys.argv[1:]; census(int(a[0]), int(a[1]), int(a[2]), int(a[3]) if len(a) > 3 else 97, pat2=(len(a) > 4 and a[4] == 'pat2'))
