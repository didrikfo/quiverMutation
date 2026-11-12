"""Round 125 (toolsmith): the opt-in library refusal `search.mutationSearchDepthFirst(..., refuseNonTilting=True)` (J != 0 steps refused,
test `mutation.isTiltingAtVertex`) measured against what it must not harm and checked against the round-033 J label.

  agree ROOT DEPTH       : guarded F/R BFS from an LNA (as 122/toolsmith_jxsearch, ISO off): at every gate-admitted step compare
                           mutation.isTiltingAtVertex with `perI` (round 033, J != 0 iff non-empty); print the 2x2 table.
  f041 N [LO HI]         : F-041 pairs at length N (pairs as 123/skeptic_f041.py), depth 3 per side, alsoDual: search.meetingPoints-style key
                           sets with the option off and on; report reached keys removed, meetings lost, and time both ways.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_refuse.py agree 003033030 3
                   timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_refuse.py f041 8"""
import sys, time, copy, collections
A = sys.argv[:]; sys.path.insert(0, '.'); sys.argv = ['x', 'none']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))   # perI, nk, mutation, procedure, reduction, ap, search, fingerprint, pathAlgebra, nx
sys.argv = A

def agree(root, D):
    r = nk.LinearNakayamaAlgebra.fromClassName(root); base = search._coxeterKeyOrNone(r)
    seen, fr, T = set(), [copy.deepcopy(r), pathAlgebra.dualPathAlgebra(r)], collections.Counter(); t0 = time.time()
    for _ in range(D):
        nf = []
        for alg in fr:
            if list(nx.simple_cycles(alg.quiver)): continue
            for kind in 'FR':
                a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
                for v in sorted(a.quiver.nodes):
                    if not mutation.mutationIsPossibleAtVertex(a, v): continue
                    J = bool(perI(a, v)); tl = mutation.isTiltingAtVertex(a, v); T[('J!=0' if J else 'J=0', 'tilting' if tl else 'not')] += 1
                    ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
                    if kind == 'R': ch = pathAlgebra.dualPathAlgebra(ch)
                    if search._coxeterKeyOrNone(ch) != base: continue
                    k = fingerprint.canonicalKey(ch)
                    if k in seen: continue
                    seen.add(k); nf.append(ch)
        fr = nf
    print('agree %s depth %d: gate-admitted steps %d, (J, isTiltingAtVertex) %s, disagreements %d, %.0fs' % (
        root, D, sum(T.values()), dict(T), T[('J!=0', 'tilting')] + T[('J=0', 'not')], time.time() - t0))

def reached(alg, depth, flag):
    found = {}
    for dual in (False, True):
        start = pathAlgebra.dualPathAlgebra(alg) if dual else copy.deepcopy(alg)
        def visit(q, path, dual=dual):
            if dual: q = pathAlgebra.dualPathAlgebra(q)
            k = search.quiverKey(q)
            if k is not None: found.setdefault(k, path[:])
        search.mutationSearchDepthFirst(start, depth, [], 'meet', printOutput=False, visitor=visit, refuseNonTilting=flag)
    return found

def f041(N, lo=0, hi=None):
    _f = open('workshop/rounds/123/skeptic_f041.py').read(); g = {'N': N}
    from quivermutation import freeMoves as fm, lnaMoves as lm
    exec(compile(_f[_f.index('def pairs():'):_f.index('def half(')], 'f041', 'exec'), dict(globals(), fm=fm, lm=lm, N=N), g)
    P = g['pairs']()[lo:hi]; tt = [0.0, 0.0]; tot = collections.Counter()
    for a, b in P:
        res = {}
        for flag in (False, True):
            t0 = time.time()
            L = reached(nk.LinearNakayamaAlgebra(N, list(a)), 3, flag); R = reached(nk.LinearNakayamaAlgebra(N, list(b)), 3, flag)
            tt[flag] += time.time() - t0; res[flag] = (set(L), set(R), len(set(L) & set(R)))
        lost = len(res[False][0] - res[True][0]) + len(res[False][1] - res[True][1])
        tot['pairs'] += 1; tot['keys off'] += len(res[False][0]) + len(res[False][1]); tot['keys removed by option'] += lost
        tot['pairs meeting off'] += res[False][2] > 0; tot['pairs meeting on'] += res[True][2] > 0
        print('  %s %s: keys %d/%d, removed %d, meetings off %d on %d' % (a, b, len(res[False][0]), len(res[False][1]), lost, res[False][2], res[True][2]), flush=True)
    print('f041 n=%d pairs %d: %s; time off %.0fs on %.0fs (x%.2f)' % (N, len(P), dict(tot), tt[0], tt[1], tt[1] / max(tt[0], 1e-9)))

if A[1] == 'agree': agree(A[2], int(A[3]))
else: f041(int(A[2]), *(int(x) for x in A[3:5]))
