"""Round 123 (skeptic), T10: the library's own merges.  F-041's meeting-in-the-middle joins (search.meetingPoints, key-guarded gate walks of
right mutations from an LNA and from its opposite, NO J filter) of an LNA and its single two-arrow deletion, at n = N, depth 3 per side.
Witness paths were stored only at n = 8 (E-156, 60 of 60 edges J = 0 by `tiltingPlus`); here they are regenerated and every edge classified:
gate, J != 0 by `perI` (the E-317 classifier's test), `tiltingPlus`, key kept, and INV = (SNF(S), signature, q_S) of every node
(skeptic_mergereplay.INV).  The B half is traversed backwards in the A -> B route; its inverse steps (forward step at v on the opposite of the
next node, as toolsmith_witness.py) are classified too.
Pairs as workshop/rounds/047/toolsmith_witness.py (`fm.derivedOrbits` with lm.ALL_MOVES; pairs the moves leave in different orbits).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py N [first last]"""
import sys, time, copy, collections
MYARGV = sys.argv[:]; N = int(MYARGV[1])
sys.path.insert(0, '.'); sys.argv = ['x', 'none']
_m = open('workshop/rounds/123/skeptic_mergereplay.py').read(); exec(compile(_m.split("PP = lambda s:")[0], 'mr', 'exec'))   # machinery + INV
from quivermutation import freeMoves as fm, lnaMoves as lm

def pairs():
    lnas, orbits = fm.derivedOrbits(N, rules=lm.ALL_MOVES, free=False, edges=True, doubles=True)
    where = {m: k for k, ms in orbits.items() for m in ms}; out = []
    for lna in lnas:
        for s, a in enumerate(lna):
            if a != 2: continue
            w = list(lna); w[s] = 0; w = tuple(w)
            if where[lna] != where[w]: out.append((lna, w))
    return out

def half(alg0, path):
    """forward edges and inverse edges of one half; rows (dir, v, gate, J, tp, kept, INV same)."""
    dual = bool(path) and path[0] < 0
    alg = pathAlgebra.dualPathAlgebra(alg0) if dual else copy.deepcopy(alg0); base = search._coxeterKeyOrNone(alg); rows = []
    for s in path:
        v = abs(s); i0 = INV(alg)
        g = mutation.mutationIsPossibleAtVertex(alg, v); J = bool(perI(alg, v)); tp = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v))
        nxt = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(alg), v))
        rows.append(('fwd', v, g, J, tp, search._coxeterKeyOrNone(nxt) in (None, base), INV(nxt) == i0))
        op = pathAlgebra.dualPathAlgebra(nxt)
        rows.append(('inv', v, mutation.mutationIsPossibleAtVertex(op, v), bool(perI(op, v)), bool(tiltingPlus(op.quiver, procedure.relationsFrom(op), v)), True, True))
        alg = nxt
    return rows, (pathAlgebra.dualPathAlgebra(alg) if dual else alg)

P = pairs(); lo, hi = (int(MYARGV[2]), int(MYARGV[3])) if len(MYARGV) > 3 else (0, len(P))
print('N', N, 'pairs', len(P), 'running', lo, hi, flush=True)
tot = collections.Counter(); t0 = time.time(); flagged = []
for i in range(lo, min(hi, len(P))):
    a, b = P[i]; X, Y = (nk.LinearNakayamaAlgebra(N, list(t)) for t in (a, b))
    ms = search.meetingPoints(X, Y, 3, alsoDual=True)
    if not ms: tot['no meeting'] += 1; print(i, ''.join(map(str, a)), ''.join(map(str, b)), 'NO MEETING', flush=True); continue
    key, pl, pr = ms[0]
    rl, el = half(X, pl); rr, er = half(Y, pr); rows = rl + rr
    tot['halves end at meeting key'] += (search.quiverKey(el) == key == search.quiverKey(er)); tot['pairs with >1 meeting'] += len(ms) > 1   # referee item 5
    fw = [r for r in rows if r[0] == 'fwd']; iv = [r for r in rows if r[0] == 'inv']
    nJ = sum(r[3] for r in fw); nJi = sum(r[3] for r in iv); ninv = sum(not r[6] for r in fw)
    tot['merges'] += 1; tot['fwd edges'] += len(fw); tot['fwd J!=0'] += nJ; tot['inv edges'] += len(iv); tot['inv J!=0'] += nJi
    tot['fwd not gate'] += sum(not r[2] for r in fw); tot['fwd tp false'] += sum(not r[4] for r in fw); tot['inv not gate'] += sum(not r[2] for r in iv)
    tot['INV changes'] += ninv; tot['merges with J!=0 (fwd or inv)'] += (nJ + nJi > 0); tot['ends INV equal'] += INV(X) == INV(Y)
    if nJ + nJi or ninv: flagged.append((i, a, b, pl, pr, nJ, nJi, ninv))
    print('%d %s -> %s halves %s | %s  fwd J!=0 %d  inv J!=0 %d  INV changes %d  %.0fs' % (i, ''.join(map(str, a)), ''.join(map(str, b)), pl, pr, nJ, nJi, ninv, time.time() - t0), flush=True)
print('TOTAL', dict(tot)); print('FLAGGED', flagged); print('time %.0fs' % (time.time() - t0))
