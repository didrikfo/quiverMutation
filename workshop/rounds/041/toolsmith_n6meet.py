"""Round 041 (toolsmith): n = 6 meet-in-the-middle, extension of rounds/038/toolsmith_n6meet.py.
 --tilting-only   BFS steps are filtered by tiltingPlus (Ladkani 2.3(c), J = 0) as well as the gate and the Coxeter key guard.
 --control        positive control: a node z with two tilting parents p1, p2 in the LNA tilting BFS; BFS(p1) and BFS(p2) must share z.
 --reverse        also search BACKWARD from the hits: BFS on the opposite algebras, forward tilting steps, keys taken of the opposite of each node
                  (reverse of a tilting step is assumed to be the forward step on the opposite algebra; checked by --revcontrol).
 --revcontrol     check on LNA-side edges A -> B: opposite(forward_tilting(opposite(B), v)) contains A.
 --plan           print sizes only (seconds per BFS = --secs, default 20) and exit
Every BFS reports closed (frontier exhausted) vs cap (time ran out).  Exit code 2 if --budget-hours is spent.
Usage: .venv/bin/python workshop/rounds/041/toolsmith_n6meet.py --tilting-only --lna-secs 150 --hit-secs 60 [--hits 0,1,..] [--control] [--reverse]"""
import sys, time, argparse
ap2 = argparse.ArgumentParser()
ap2.add_argument('--tilting-only', action='store_true'); ap2.add_argument('--control', action='store_true')
ap2.add_argument('--reverse', action='store_true'); ap2.add_argument('--revcontrol', action='store_true'); ap2.add_argument('--plan', action='store_true')
ap2.add_argument('--lna-secs', type=float, default=150); ap2.add_argument('--hit-secs', type=float, default=30)
ap2.add_argument('--hits', default='all'); ap2.add_argument('--budget-hours', type=float, default=0.15)
A2 = ap2.parse_args()
T0 = time.time(); BUDGET = A2.budget_hours * 3600
if A2.plan: A2.lna_secs = A2.hit_secs = 20
src = open('workshop/rounds/038/toolsmith_n6close.py').read().split("budget = args")[0]
sys.argv = ['x', '--plan', '1']; exec(compile(src, 'c', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]
exec(compile(h15, 'h015', 'exec'))
key0 = order[0]
def spent():
    if time.time() - T0 > BUDGET: print('BUDGET SPENT'); sys.exit(2)
def bfs(starts, secs, tilting, opp=False):
    """seen: key -> (algebra, [(parentkey, v)], depth). With opp, nodes are opposite algebras and keys are canonicalKey of the opposite (true) algebra."""
    seen = {}; fr = []; t0 = time.time(); stats = Counter(); capped = False
    def tk(x): return fingerprint.canonicalKey(pathAlgebra.dualPathAlgebra(x) if opp else x)
    def add(x, par, v, d):
        ck = tk(x)
        if ck is None: stats['nokey'] += 1; return
        if ck in seen:
            if par is not None and (par, v) not in seen[ck][1]: seen[ck][1].append((par, v))
            return
        seen[ck] = (x, [] if par is None else [(par, v)], d); fr.append((x, ck, d))
    for x in starts: add(x, None, None, 0)
    while fr and not capped:
        cur, fr[:] = fr[:], []
        for j, (alg, ck0, d) in enumerate(cur):
            if time.time() - t0 > secs: capped = True; break
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                stats['gate'] += 1
                if tilting and not tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v): stats['J!=0 dropped'] += 1; continue
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == key0: add(ch, ck0, v, d + 1)
    return seen, not capped, stats, len(fr)
def rep(label, S, closed, stats, left):
    print('%s: seen %d, %s%s, gate-admitted %d, J!=0 dropped %d' % (label, len(S), 'CLOSED' if closed else 'CAP HIT (not closed)',
          '' if closed else ', frontier left %d' % left, stats['gate'], stats['J!=0 dropped']), flush=True)
til = A2.tilting_only
print('mode: tilting-only=%s' % til, flush=True)
S1, c1, s1, l1 = bfs(classes[key0], A2.lna_secs, til); rep('LNA side', S1, c1, s1, l1); spent()
if A2.control and til:
    z = None
    for k, (x, ps, d) in sorted(S1.items(), key=lambda kv: kv[1][2]):
        pk = {p for p, v in ps}
        if len(pk) >= 2 and all(S1[p][2] >= 1 for p in pk): z = k; break
    if z is None: print('CONTROL: no node with two distinct tilting parents at depth>=1 found'); 
    else:
        p1, p2 = sorted({p for p, v in S1[z][1]}, key=lambda p: S1[p][2])[:2]
        X, c, s, l = bfs([S1[p1][0]], A2.lna_secs / 2, True); Y, c2, s2, l2 = bfs([S1[p2][0]], A2.lna_secs / 2, True)
        print('CONTROL z depth %d (LNA-side), parents at depths %d, %d; p2 in BFS(p1): %s' % (S1[z][2], S1[p1][2], S1[p2][2], p2 in X))
        rep('  BFS(p1)', X, c, s, l); rep('  BFS(p2)', Y, c2, s2, l2)
        print('CONTROL RESULT: z in both = %s; shared = %d' % (z in X and z in Y, len(set(X) & set(Y))), flush=True)
if A2.revcontrol and til:
    ok = bad = 0; t1 = time.time()
    for k, (x, ps, d) in list(S1.items()):
        if d == 0 or time.time() - t1 > A2.lna_secs / 2: continue
        for p, v in ps[:1]:
            B = x; Bop = pathAlgebra.dualPathAlgebra(B)
            found = False
            if not list(nx.simple_cycles(Bop.quiver)):
                for w in sorted(Bop.vertices()):
                    if not mutation.mutationIsPossibleAtVertex(Bop, w) or not tiltingPlus(Bop.quiver, procedure.relationsFrom(Bop), w): continue
                    raw = mutation.quiverMutationAtVertex(Bop, w)
                    if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                    ch = reduction.reducePathAlgebra(raw)
                    if fingerprint.canonicalKey(pathAlgebra.dualPathAlgebra(ch)) == p: found = True; break
            ok += found; bad += (not found)
    print('REVCONTROL: LNA-side edges A->B tested %d, A recovered as opposite(step(opposite(B))) in %d, not in %d' % (ok + bad, ok, bad), flush=True)
sel = range(len(hits)) if A2.hits == 'all' else [int(x) for x in A2.hits.split(',')]
starts = {j: build([tuple(x) for x in hits[j][0]], hits[j][1], list(range(1, 7))) for j in sel}
inl = sum(1 for k in hitkeys if k in S1); print('hits inside the LNA-side set:', inl, flush=True)
tot = Counter(); closedhits = 0
for j in sel:
    spent()
    S2, c2, s2, l2 = bfs([starts[j]], A2.hit_secs, til); m = set(S1) & set(S2); closedhits += c2
    rep('hit %d (v=%d) forward' % (j, hits[j][3]), S2, c2, s2, l2); print('   shared with LNA side: %d' % len(m), flush=True)
    tot['fwd'] += len(m)
    if A2.reverse:
        R, cr, sr, lr = bfs([pathAlgebra.dualPathAlgebra(starts[j])], A2.hit_secs, til, opp=True); mr = set(S1) & set(R)
        rep('hit %d (v=%d) REVERSE(J=0 only)' % (j, hits[j][3]), R, cr, sr, lr); print('   shared with LNA side: %d' % len(mr), flush=True)
        tot['rev'] += len(mr)
print('SUMMARY tilting-only=%s hits=%d forward-closed hits=%d shared fwd=%d rev=%d LNA closed=%s' % (til, len(sel), closedhits, tot['fwd'], tot['rev'], c1))
