"""Round 097 skeptic, T10/T5: do E-235 / E-248-E-249 cover every step the REAL search takes from acyclic parents that carry a
degree-1 relation arrow omega: i ~> t parallel to an arrow i -> t (the pattern that broke E-248)?
Usage (repo root):  .venv/bin/python workshop/rounds/097/skeptic_omega.py NMIN NMAX PARENTS_PER_N DEPTH [seed]
Per parent: the real search (search.mutationSearchDepthFirst, default coxeterGuard) with a visitor; every visited node with
mutationVertices s is the child of the node s[:-1] at vertex s[-1].  Per step (parent acyclic, because only acyclic nodes expand):
  H   : E-235 hypothesis 'no loop at i' = parent Gabriel quiver has no oriented cycle (checked, trivially true when acyclic)
  PAT : pattern of E-248 = a parent relation with source i and target t and an arrow i -> t exist (omega parallel to alpha)
  MIN : parent relation set minimal? (no relation in the ideal of the others; a redundant one needs a degree-2 generator)
  child cyclic? ; and, for an acyclic child, the independent End(T) test: dim Hom_K(T_a,T_b) (round-050 code, T = F tilting complex)
        against the labelled Cartan matrix of the repo child (either orientation, as E-181)."""
import sys, random, io, contextlib, itertools
sys.path.insert(0, '.')
from collections import Counter
import networkx as nx
_st = open('workshop/rounds/050/skeptic_tilt.py').read()
exec(_st[:_st.index('def stepTest')] + _st[_st.index('def stepTest'):])
from quivermutation import search, procedure, arrowPaths as ap

def gen_parent(rng, n, want_pat=True, pat2=False):
    q = nx.MultiDiGraph(); q.add_nodes_from(range(n))
    for a in range(n):
        for b in range(a + 1, n):
            if rng.random() < 0.45: q.add_edge(a, b)
    paths = {}
    rels = []
    # plant: i<j<t, arrows i->j, j->t, and arrow i->t, zero relation on the 2-path (omega || alpha)
    chain3 = None
    if pat2:
        i, j, k, t = sorted(rng.sample(range(n), 4)); want_pat = False
        for (x, y) in ((i, j), (j, k), (k, t), (i, t)):
            if not q.has_edge(x, y): q.add_edge(x, y)
        chain3 = (i, j, k, t)
    if want_pat:
        i, j, t = sorted(rng.sample(range(n), 3))
        for (x, y) in ((i, j), (j, t), (i, t)):
            if not q.has_edge(x, y): q.add_edge(x, y)
    # extra zero relations of length 2 (random), commutativity of two parallel 2-paths
    two = []
    for b in q.nodes:
        for e1 in list(q.in_edges(b, keys=True)):
            for e2 in list(q.out_edges(b, keys=True)):
                two.append((e1, e2))
    rng.shuffle(two)
    chosen = []
    if chain3:
        i, j, k, t = chain3
        for (x, y, z) in ((i, j, k), (j, k, t)):
            chosen.append(((x, y, min(q[x][y])), (y, z, min(q[y][z]))))
    if want_pat:
        # make sure the planted pair is first
        cand = [(e1, e2) for (e1, e2) in two if e1[0] == i and e2[1] == t and e1[1] == j]
        if cand: chosen.append(cand[0])
    for pair in two:
        if len(chosen) >= 1 + rng.randint(0, 2): break
        if pair not in chosen and rng.random() < 0.5: chosen.append(pair)
    for (e1, e2) in chosen:
        rels.append(ap.combination([((e1, e2), 1)]))
    # sometimes a commutativity relation between two 2-paths with the same ends
    ends = {}
    for (e1, e2) in two: ends.setdefault((e1[0], e2[1]), []).append((e1, e2))
    for k, v in ends.items():
        if len(v) >= 2 and rng.random() < 0.4 and all(p not in chosen for p in v[:2]):
            rels.append(ap.combination([(v[0], 1), (v[1], -1)]))
    # sometimes zero relations of arrow length 3 (a 2-path then an arrow), so overlaps / Ext^3 occur
    for (e1, e2) in two:
        for e3 in q.out_edges(e2[1], keys=True):
            if rng.random() < 0.08: rels.append(ap.combination([((e1, e2, e3), 1)]))
    if not rels: return None
    if minimal(q, rels) is not True: return None      # roots are minimal presentations; descendants are reported as found
    return procedure.toPathAlgebra(q, rels)

def acyc(q): return not ap.hasOrientedCycle(q)

def minimal(q, rels):
    for k, r in enumerate(rels):
        others = rels[:k] + rels[k + 1:]
        try:
            if ap.isInIdeal(q, others, r): return False
        except Exception: return None
    return True

def run(nmin, nmax, per, depth, seed):
    rng = random.Random(seed); cnt = Counter(); bad = []
    for n in range(nmin, nmax + 1):
        done = 0; tries = 0
        while done < per and tries < 20 * per:
            tries += 1
            P0 = gen_parent(rng, n)
            if P0 is None: continue
            seen = {}
            def visitor(pa, seq):
                import copy
                seen.setdefault(tuple(seq), copy.deepcopy(pa))
            with contextlib.redirect_stdout(io.StringIO()):
                search.mutationSearchDepthFirst(P0, depth, printOutput=False, visitor=visitor)
            done += 1; cnt['n=%d parents' % n] += 1
            for seq, ch in seen.items():
                if not seq: continue
                par = seen.get(seq[:-1]); v = seq[-1]
                if par is None: cnt['parent not recorded'] += 1; continue
                q = par.quiver; rels = procedure.relationsFrom(par)
                H = acyc(q)
                pat = any(ap.source(r) == v and any(a[1] == ap.target(r) for a in ap.arrowsOutOf(q, v)) for r in rels)
                mn = minimal(q, rels)
                chq = ch.quiver; cyc = not acyc(chq)
                key = ('n=%d' % n, 'H' if H else 'noH', 'PAT' if pat else 'nopat', 'min' if mn else ('NONMIN' if mn is False else 'min?'), 'cyclic child' if cyc else 'acyclic child')
                cnt[key] += 1
                if not H: bad.append(('noH', seq))
                if cyc: continue
                try:
                    Vs, res = stepTest(par, v, -1)
                    C = ap.cartanMatrix(chq, procedure.relationsFrom(ch), exact=True)
                except Exception as e:
                    cnt[key + ('ERR ' + type(e).__name__,)] += 1; continue
                Hm = res[0]; N = len(Vs)
                eq = all(C[a][b] == Hm[a][b] for a in range(N) for b in range(N)) or all(C[a][b] == Hm[b][a] for a in range(N) for b in range(N))
                J = sum(map(sum, res[-1])) + sum(map(sum, res[1]))
                cnt[key + ('End=child' if eq else 'MISMATCH', 'J=0' if J == 0 else 'J!=0')] += 1
                if not eq: bad.append(('mismatch', n, seq, [list(x) for x in Hm], C))
    for k in sorted(cnt, key=str): print(cnt[k], k)
    print('bad:', len(bad)); [print(b) for b in bad[:5]]

if __name__ == '__main__':
    a = sys.argv[1:]
    run(int(a[0]), int(a[1]), int(a[2]), int(a[3]), int(a[4]) if len(a) > 4 else 97)
