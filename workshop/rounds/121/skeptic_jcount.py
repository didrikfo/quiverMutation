"""Round 121 (skeptic): T10 / H-015 power pre-check at n = 9 (F-010 key group).
mode `ball ROOT DEPTH`: rebuild E-306's ball (E-212 `moves`, J = 0 key-kept tiltingPlus) and, at every node of depth <= DEPTH, tally EVERY
  gate-admissible perI (J != 0) F/R candidate by outcome (kept / moved / illegal / verts / err) -- E-306's script dropped a kept candidate
  silently if `child` raised; here such a drop is counted -- and the q_S invariant (signature, |det S|, q_S histogram of S = C + C^T) of each
  key-kept child, labelled A (= 3060000 half), B (= 3304000 half) or OTHER.
mode `replay`: rebuild the parents of all 19 recorded key-kept perI J != 0 steps of the F-010 balls (E-233, E-234 x10, E-239, E-244 x6, E-254) from
  their recorded paths, check every path move is a `moves` move, check the step is perI-flagged and key-kept, and label the child's q_S.
  Planted control: each step's child compared against the OTHER half's root invariant (must say NO).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/121/skeptic_jcount.py ball 3060000 4 | replay"""
import sys, os, time, collections
MODE = sys.argv[1:]
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
_sj = open('workshop/rounds/082/skeptic_join.py').read(); exec(_sj[_sj.index('def child('):_sj.index('def struct(')])
sys.path.insert(0, "workshop/rounds/117")
from fractions import Fraction
import numpy as np, sympy as sp
from quivermutation import invariants
_q = open('workshop/rounds/117/toolsmith_qsball.py').read(); exec(_q[_q.index('_cache = {}'):_q.index('def ballrun(')])   # qinv, unchanged
INV = {r: qinv(nk.LinearNakayamaAlgebra.fromClassName(r)) for r in ('3060000', '3304000', '2400230')}
assert INV['3304000'] == INV['2400230'] and INV['3060000'] != INV['3304000']
def half(alg):
    q = qinv(alg); return 'A' if q == INV['3060000'] else 'B' if q == INV['3304000'] else 'OTHER'

def jscan(alg):
    """every gate-admissible perI F/R candidate: (kind, v, tag, child or None)."""
    V = sorted(alg.quiver.nodes); out = []
    for kind in 'FR':
        a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(a, v) or not perI(a, v): continue
            try:
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): out.append((kind, v, 'illegal', None)); continue
                c = reduction.reducePathAlgebra(raw)
                if sorted(c.vertices()) != V: out.append((kind, v, 'verts', None)); continue
                if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
                k = search._coxeterKeyOrNone(c)
                out.append((kind, v, 'keyNone' if k is None else ('kept' if k == base else 'moved'), c))
            except Exception as e: out.append((kind, v, 'err ' + type(e).__name__, None))
    return out

if MODE[0] == 'ball':
    root, D = MODE[1], int(MODE[2]); t0 = time.time()
    r = nk.LinearNakayamaAlgebra.fromClassName(root); seen = {fingerprint.canonicalKey(r): r}; front = [r]
    tags = collections.Counter(); kept = []
    for d in range(D + 1):                    # J scan at depths 0..D (E-306's LASTJ=1 default), J = 0 expansion to depth D
        for i, x in enumerate(front):
            for kind, v, tag, c in jscan(x):
                tags[(d, tag)] += 1
                if tag == 'kept':
                    try: c2 = child(x, kind, v); built = True
                    except Exception as e: c2 = c; built = 'child() raised ' + type(e).__name__
                    kept.append((d, i, kind, v, built, half(x), half(c2), qinv(x) == qinv(c2)))
        if d == D: break
        nxt = []
        for x in front:
            for kind, v, c in moves(x):
                k = fingerprint.canonicalKey(c)
                if k not in seen: seen[k] = c; nxt.append(c)
        front = nxt
        print(root, 'depth', d + 1, 'nodes', len(seen), 'front', len(front), '%.0fs' % (time.time() - t0), flush=True)
    print('ROOT', root, 'DEPTH', D, 'nodes', len(seen), 'root half', half(r))
    print('perI candidate outcomes by parent depth:', dict(sorted(tags.items())))
    print('KEY-KEPT J != 0 steps:', len(kept))
    for k in kept: print('  parent depth %d idx %d step %s%d built=%s parent half %s child half %s q_S equal %s' % k)
    print('time %.0fs' % (time.time() - t0))

elif MODE[0] == 'replay':
    PP = lambda s: [(t[0], int(t[1:])) for t in s.split()]
    STEPS = [  # (record, root, path to parent, step)
        ('E-233', '3060000', 'F1 R9 R8 R8', 'R8'),
        ('E-234 1275', '3060000', 'F1 F4 R9 R8 R8', 'R8'), ('E-234 1433', '3060000', 'F1 F5 R9 R8 R8', 'R8'),
        ('E-234 1746', '3060000', 'F1 R5 R9 R8 R8', 'R8'), ('E-234 1818', '3060000', 'F1 R6 R9 R8 R8', 'R8'),
        ('E-234 1915', '3060000', 'F1 R8 R9 F8 R9', 'R9'), ('E-234 1932', '3060000', 'F1 R8 R9 R9 R8', 'R8'),
        ('E-234 1970', '3060000', 'F1 R9 F1 R8 R8', 'R8'), ('E-234 2003', '3060000', 'F1 R9 F3 R8 R8', 'R8'),
        ('E-234 4551', '3060000', 'R8 R9 R9 F1 R8', 'R8'), ('E-234 4722', '3060000', 'R9 R8 R8 F1 F8', 'R1'),
        ('E-239 2026', '3304000', 'F2 R8 R9 R7 R7', 'R7'),
        ('E-244 4606', '3304000', 'F2 F1 R8 R9 R7 R7', 'R7'), ('E-244 8063', '3304000', 'F2 R8 R7 R9 F7 R9', 'R9'),
        ('E-244 8170', '3304000', 'F2 R8 R9 F2 R7 R7', 'R7'), ('E-244 8196', '3304000', 'F2 R8 R9 F4 R7 R7', 'R7'),
        ('E-244 8427', '3304000', 'F2 R9 R8 R7 R9 R7', 'R7'), ('E-244 11123', '3304000', 'R8 R9 R7 R7 F2 F7', 'R2'),
        ('E-254', '2400230', 'F1 F1 F2 F3 R9 R7 R9 F3', 'F3')]   # added round 121 response (referee: 19th step, 2400230 is half B)
    t0 = time.time(); tally = collections.Counter()
    for rec, root, path, step in STEPS:
        x = nk.LinearNakayamaAlgebra.fromClassName(root); ok = True
        for kind, v in PP(path):
            mv = [c for k_, v_, c in moves(x) if (k_, v_) == (kind, v)]
            if not mv: ok = False; break
            x = mv[0]
        (kind, v), = PP(step)
        js = [t for t in jscan(x) if (t[0], t[1]) == (kind, v)] if ok else []
        tag = js[0][2] if js else 'not perI-flagged'
        c = child(x, kind, v) if ok else None
        hp, hc = (half(x), half(c)) if ok else ('-', '-')
        other = INV['3304000'] if root == '3060000' else INV['3060000']
        pm = PP(path); shape = ('step=last move' if pm[-1] == (kind, v) else 'step in path' if (kind, v) in pm else 'step not in path')
        tally[('shape', shape)] += 1   # printed shape key: does the step repeat a move (same kind, same vertex) of the path?
        planted_no = (qinv(c) != other) if ok else None       # planted cross step: child vs the other half's root; must be True (NO)
        res = 'YES (q_S equal)' if ok and qinv(c) == qinv(x) else 'NO' if ok else 'path broken'
        tally[(tag, res, planted_no)] += 1
        print('%-12s root %s path %-24s step %s  shape %-16s path moves ok %s  tag %-8s parent %s child %s  q_S test: %s  planted cross step NO: %s'
              % (rec, root, path, step, shape, ok, tag, hp, hc, res, planted_no), flush=True)
    print('TALLY', dict(tally), 'time %.0fs' % (time.time() - t0))
