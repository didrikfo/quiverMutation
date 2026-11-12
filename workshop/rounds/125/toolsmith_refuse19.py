"""Round 125 (toolsmith), answer to the referee's item 1: the 19 audited n = 9 F-010 J != 0 steps (E-316) under refuseNonTilting.

For each recorded join (the 14 printed in 123/skeptic_mergereplay.py `records` + its 5 `extra` = 19 steps):
  - replay the parent path, the step, the child path and the other side's root path with the round-123 classifier (gate, J, tiltingPlus, key);
  - check `mutation.isTiltingAtVertex` at the audited step (False = the option refuses it);
  - check every OTHER step of the join is tilting (so the option refuses only the audited step);
  - option-on detour: guarded F/R BFS from the parent P using only steps the option admits (gate, key kept, isTiltingAtVertex), up to
    depth DMAX; report the least depth at which a node isomorphic to the child X appears (E-316: X is in P's J = 0 component).
    The join then survives with the option on, its audited step replaced by that detour (length recorded - 1 + detour).
usage (repo root): DMAX=3 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_refuse19.py"""
import sys, os, time, copy
sys.path.insert(0, '.'); A = sys.argv[:]; sys.argv = ['x', 'none']
_m = open('workshop/rounds/123/skeptic_mergereplay.py').read()
exec(compile(_m[_m.index('import sys, time'):_m.index("if MODE == 'records':")].replace('MODE = sys.argv[1]; ', ''), 'mr', 'exec'))
sys.argv = A
RECS = eval(_m[_m.index('RECS = [') + 7:_m.index("]\n    out = []; kids") + 1].replace('A)', "'F1 F7 F6 F5 F4')"))
EXTRA = eval(_m[_m.index('EXTRA = [') + 8:_m.index("]\n    for rec, pp, st, cp, rp") + 1])
RECS += [(r, '3060000', pp, st, cp, '3060000', rp) for r, pp, st, cp, rp in EXTRA]
# E-255 8427: child iso to E-244 8063's child (skeptic_mergereplay `records`), so its join is 8063's child path and root path
RECS += [('E-255 8427', '3304000', 'F2 R9 R8 R7 R9 R7', 'R7', 'R2 R4 R5', '3304000', 'F2 F1 F6')]
# 'E-233 X~P' re-uses E-233's step (second join of the same step): 20 joins, 19 distinct audited steps
DMAX = int(os.environ.get('DMAX', '3'))

def detour(P, X, pk0):
    sx = sig(X); seen = {fingerprint.canonicalKey(P)}; fr = [P]
    for d in range(1, DMAX + 1):
        nf = []
        for alg in fr:
            for kind in 'FR':
                a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
                for v in sorted(a.quiver.nodes):
                    if not mutation.mutationIsPossibleAtVertex(a, v) or not mutation.isTiltingAtVertex(a, v): continue
                    c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
                    if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
                    if search._coxeterKeyOrNone(c) != pk0: continue
                    k = fingerprint.canonicalKey(c)
                    if k in seen: continue
                    seen.add(k)
                    if sig(c) == sx and iso(c, X): return d, len(seen)
                    nf.append(c)
        fr = nf
    return None, len(seen)

t0 = time.time(); pk0 = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
tot = dict(refused=0, others_tilting=0, detour=0)
print('%-12s %-6s %-9s %-14s %-8s %s' % ('record', 'step', 'tilting?', 'other steps', 'detour', 'nodes'))
for rec, root, pp, st, cp, rr, rp in RECS:
    rows = []; r = nk.LinearNakayamaAlgebra.fromClassName(root)
    P = walk(r, pp, pk0, 'parent path', rows)
    (kind, v), = PP(st); a = P if kind == 'F' else pathAlgebra.dualPathAlgebra(P)
    tl = mutation.isTiltingAtVertex(a, v)
    X = walk(P, st, pk0, 'STEP', rows); walk(X, cp, pk0, 'child path', rows)
    R2 = nk.LinearNakayamaAlgebra.fromClassName(rr) if rr != 'PARENT' else P; walk(R2, rp, pk0, 'root path', rows)
    others = [t for t in rows if t[0] != 'STEP']
    oth_ok = all(t[2] == 'J=0' and t[3] == 'tp' for t in others)
    d, nn = detour(P, X, pk0)
    tot['refused'] += (not tl); tot['others_tilting'] += oth_ok; tot['detour'] += d is not None
    print('%-12s %-6s %-9s %-14s %-8s %d' % (rec, st, tl, 'all J=0 tp' if oth_ok else 'SOME J!=0', d, nn), flush=True)
print('of %d audited steps: refused by the option %d; joins whose other steps are all J = 0 tilting %d; child reached from parent by an '
      'option-admitted detour of depth <= %d: %d; %.0fs' % (len(RECS), tot['refused'], tot['others_tilting'], DMAX, tot['detour'], time.time() - t0))
