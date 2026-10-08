"""Round 043 (toolsmith): reverse positive control at depth 2 and diagnosis of the ~12% lost reverse edges (E-142).
Uses the prelude of rounds/041/toolsmith_n6meet.py (LNA class 0 at n = 6, gate, Coxeter key guard, tiltingPlus).
 Part A (--diag): collect forward tilting edges A -(v)-> B on the LNA side; for each, look for a vertex w of opposite(B) whose
          UNFILTERED mutation has opposite == A (up to canonical key), then say which filter (gate / J / illegal relation / key guard / none)
          rejects it. Also records whether w == v.
 Part B (--control): a two-step forward path A -> B -> C (A a start LNA, C at depth 2); depth-limited reverse search from C
          (opposite algebras, forward tilting steps, keys of the opposite) must reach A within 2 steps; also a depth-1 control and a negative
          (reverse from an unrelated node must not reach A).
Exit code 2 if --budget-hours spent.  Usage: .venv/bin/python workshop/rounds/043/toolsmith_revcontrol.py --diag --control --secs 60"""
import sys, time, argparse
ap2 = argparse.ArgumentParser()
ap2.add_argument('--diag', action='store_true'); ap2.add_argument('--control', action='store_true')
ap2.add_argument('--secs', type=float, default=60); ap2.add_argument('--budget-hours', type=float, default=0.15)
ap2.add_argument('--depth', type=int, default=2); ap2.add_argument('--maxedges', type=int, default=3000); ap2.add_argument('--npaths', type=int, default=10)
A2 = ap2.parse_args()
T0 = time.time(); BUDGET = A2.budget_hours * 3600
src = open('workshop/rounds/038/toolsmith_n6close.py').read().split("budget = args")[0]
sys.argv = ['x', '--plan', '1']; exec(compile(src, 'c', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]
exec(compile(h15, 'h015', 'exec'))
key0 = order[0]
def spent():
    if time.time() - T0 > BUDGET: print('BUDGET SPENT'); sys.exit(2)
dual = pathAlgebra.dualPathAlgebra
def ck(x, opp=False): return fingerprint.canonicalKey(dual(x) if opp else x)
def step(alg, v, tilting=True, guard=True):
    """returns (child or None, reason)"""
    if not mutation.mutationIsPossibleAtVertex(alg, v): return None, 'gate'
    if tilting and not tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v): return None, 'J'
    raw = mutation.quiverMutationAtVertex(alg, v)
    if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): return None, 'illegal'
    ch = reduction.reducePathAlgebra(raw)
    if guard and search._coxeterKeyOrNone(ch) != key0: return None, 'keyguard'
    return ch, 'ok'
def bfs(starts, secs, maxd, opp=False):
    seen = {}; fr = []; t0 = time.time(); stats = Counter()
    def add(x, par, v, d):
        k = ck(x, opp)
        if k is None: stats['nokey'] += 1; return
        if k in seen:
            if par is not None and (par, v) not in seen[k][1]: seen[k][1].append((par, v))
            return
        seen[k] = (x, [] if par is None else [(par, v)], d); fr.append((x, k, d))
    for x in starts: add(x, None, None, 0)
    capped = False
    while fr and not capped:
        cur, fr[:] = fr[:], []
        for alg, k0, d in cur:
            if d >= maxd: continue
            if time.time() - t0 > secs: capped = True; break
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                ch, why = step(alg, v)
                stats[why] += 1
                if ch is not None: add(ch, k0, v, d + 1)
    return seen, (not capped), stats
S1, c1, s1 = bfs(classes[key0], A2.secs, 99)
print('LNA side (tilting-only): seen %d closed=%s; depth hist %s' % (len(S1), c1, sorted(Counter(d for _, _, d in S1.values()).items())), flush=True)
spent()
if A2.diag:
    cat = Counter(); sameW = Counter(); n = 0
    for k, (B, ps, d) in S1.items():
        if d == 0: continue
        for p, v in ps[:1]:
            if n >= A2.maxedges: break
            n += 1; Bop = dual(B)
            if list(nx.simple_cycles(Bop.quiver)): cat['opp has cycle'] += 1; continue
            full = None; allk = {}
            for w in sorted(Bop.vertices()):
                # unfiltered mutation (gate only needed to build)
                if not mutation.mutationIsPossibleAtVertex(Bop, w): continue
                try:
                    raw = mutation.quiverMutationAtVertex(Bop, w); ch = reduction.reducePathAlgebra(raw)
                    kk = ck(ch, True)
                except Exception as e: continue
                if kk == p: full = (w, raw, ch); break
            if full is None:
                cat['no w gives A (A is not an opposite-mutation of B at all)'] += 1
                # what does the same-vertex opposite step give instead?
                if mutation.mutationIsPossibleAtVertex(Bop, v):
                    A1 = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(Bop, v)); k1 = ck(A1, True)
                    cat['   same-vertex opposite step gives A2 != A; A2 has class key: %s; A2 in LNA BFS: %s' % (search._coxeterKeyOrNone(A1) == key0, k1 in S1)] += 1
                    a1 = len(list(A1.quiver.edges())); a0 = len(list(S1[p][0].quiver.edges()))
                    cat['   arrows(A2) - arrows(A) = %+d' % (a1 - a0)] += 1
                continue
            w, raw, ch = full
            sameW['w == v' if w == v else 'w != v'] += 1
            if not mutation.mutationIsPossibleAtVertex(Bop, w): why = 'gate'
            elif not tiltingPlus(Bop.quiver, procedure.relationsFrom(Bop), w): why = 'J'
            elif any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): why = 'illegal'
            elif search._coxeterKeyOrNone(ch) != key0: why = 'keyguard'
            else: why = 'ok (recovered)'
            cat[why] += 1
    print('DIAG edges tested %d' % n); [print('  %-60s %d' % (a, b)) for a, b in cat.most_common()]
    print('  vertex relation:', dict(sameW), flush=True)
    spent()
if A2.control:
    # two-step forward paths A -> B -> C, A depth 0
    paths = []
    def back(k, d):  # all backward chains of forward edges from k at depth d down to a depth-0 node
        if d == 0: return [(k,)]
        return [c + (k,) for pk, v in S1[k][1] if S1[pk][2] == d - 1 for c in back(pk, d - 1)]
    for kc, (C, psc, dc) in S1.items():
        if dc != A2.depth: continue
        paths += [(c[0], c[1], c[-1]) for c in back(kc, dc)[:3]]
        if len(paths) >= 600: break
    seenp = set(); res = []; import random; random.seed(0)
    paths = [p for p in paths if not (p in seenp or seenp.add(p))]; random.shuffle(paths)
    print('two-step forward paths available (sample of C at depth 2): %d' % len(paths), flush=True)
    for ka, kb, kc in paths[:A2.npaths]:
        spent()
        R, cl, st = bfs([dual(S1[kc][0])], 120, A2.depth, opp=True)
        R1, _, _ = bfs([dual(S1[kc][0])], 120, A2.depth - 1, opp=True)
        res.append((ka in R, R[ka][2] if ka in R else None, kb in R, ka in R1, len(R), len(R1), cl))
        print('path A->B->C: reverse depth<=%d from C reaches A: %s (depth %s), reaches B(first forward step): %s; depth<=%d reaches A: %s; |R2|=%d |R1|=%d R2 closed=%s' % ((A2.depth,) + res[-1][:3] + (A2.depth - 1,) + res[-1][3:]), flush=True)
    # negative control: reverse from a node not forward-reachable ... use a node in a different class key (hit) and test A
    print('CONTROL SUMMARY: A reached in %d of %d paths at depth<=%d' % (sum(r[0] for r in res), len(res), A2.depth))
