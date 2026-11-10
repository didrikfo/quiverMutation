"""Round 123 (skeptic), T10 / H-015 / H-032: replay the RECORDED join ("merge") paths of the n = 9 F-010 key group and classify every step.

A recorded join is   root --parent path--> P --step--> X --child path--> M_X  ~iso~  M_R <--root path-- root'   (E-233, E-234, E-239, E-244, E-254).
Every step (parent path, the audited step, child path, root path) is re-taken and classified exactly as the E-317 classifier
(skeptic_jtilt.py = maverick_jtilt.py) does: gate (`mutationIsPossibleAtVertex`), illegal relation, J != 0 (`perI` non-empty), `tiltingPlus`,
Coxeter key kept.  Added: the derived invariant INV = (SNF(S), signature of S, q_S histogram when det S != 0), S = C + C^T, at every node,
and the isomorphism test at the meet (skeptic_join.py `iso`).  A merge "contains a J != 0 step" if any of its steps has perI non-empty.
INV must be constant along any chain of derived equivalences; a change marks a step that is not one (E-318).

modes (run from the repository root):
  records   all recorded n = 9 joins with printed paths (13) + E-233's child-to-parent join + E-255's iso claim (8427 child ~ 8063 child)
  f037      F-037's five recorded F-only paths 34504030 -> 50505000 (n = 10, the library walk's merge)
  control   E-318's 5-step walk at n = 11 (003033030 -F3 R8 F3 F3 F6-> 030330030): must be flagged (J != 0 at step 4, INV changes)
  mkpkl DIR  write parent pickles of the five E-234 joins not printed in the record, for skeptic_join.py (regeneration; see submission)
  extra     replay the five regenerated E-234 joins (paths as printed by skeptic_join.py meet, copied into EXTRA below)
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py records|control|mkpkl DIR|extra"""
import sys, time, copy, collections
from fractions import Fraction
MODE = sys.argv[1]; ARGV = sys.argv[:]
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]; exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
_sj = open('workshop/rounds/082/skeptic_join.py').read(); exec(_sj[_sj.index('def struct('):_sj.index('def mkroot(')])   # struct, iso, sig
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import invariants

_cache = {}
def INV(alg):
    C = np.array(invariants.integerCartanMatrix(alg)).astype(int); k = (C.tobytes(), C.shape)
    if k in _cache: return _cache[k]
    n = C.shape[0]; S = sp.Matrix((C + C.T).tolist())
    snf = smith_normal_form(S, domain=sp.ZZ); dg = tuple(sorted(abs(int(snf[i, i])) for i in range(n)))
    ev = np.linalg.eigvalsh(np.array(S.tolist(), dtype=float)); sg = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
    d = int(S.det()); hh = None
    if d != 0:   # q_S histogram on the discriminant group (as toolsmith_qsball.qinv, round 117)
        Si = S.inv(); cols = [[Fraction(int(Si[i, j].p), int(Si[i, j].q)) for i in range(n)] for j in range(n)]
        kf = lambda v: tuple(x % 1 for x in v); seen = {kf([Fraction(0)] * n)}; fr = list(seen)
        while fr:
            nf = []
            for g in fr:
                for c in cols:
                    h = kf([a + b for a, b in zip(g, c)])
                    if h not in seen: seen.add(h); nf.append(h)
            fr = nf
        Sl = [[int(S[i, j]) for j in range(n)] for i in range(n)]
        hh = tuple(sorted(collections.Counter(sum(y[i] * Sl[i][j] * y[j] for i in range(n) for j in range(n)) % 2 for y in seen).items()))
    _cache[k] = (dg, sg, hh); return _cache[k]
short = lambda inv: 'SNF tail %s sig %s q %s' % (inv[0][-3:], inv[1], 'n/a(det 0)' if inv[2] is None else len(inv[2]))

def take(x, kind, v, pk0):
    """one step, classified as skeptic_jtilt.py does; returns (child or None, label dict)."""
    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    lab = dict(gate=bool(mutation.mutationIsPossibleAtVertex(a, v)))
    if not lab['gate']: return None, lab
    raw = mutation.quiverMutationAtVertex(a, v)
    lab['illegal'] = any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw))
    lab['J'] = bool(perI(a, v)); lab['tp'] = bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v))
    c = reduction.reducePathAlgebra(raw)
    if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
    lab['kept'] = search._coxeterKeyOrNone(c) == pk0
    return c, lab

PP = lambda s: [(t[0], int(t[1:])) for t in s.split()]
def walk(start, path, pk0, tag, rows):
    """replay path from start; append one row per step; return end node (None if a step is refused)."""
    x = start
    for kind, v in PP(path):
        i0 = INV(x); c, lab = take(x, kind, v, pk0)
        if c is None: rows.append((tag, kind + str(v), 'GATE REFUSES', '', '')); return None
        i1 = INV(c)
        rows.append((tag, kind + str(v), 'J!=0' if lab['J'] else 'J=0', 'tp' if lab['tp'] else 'NOT-tp',
                     ('kept' if lab['kept'] else 'MOVED') + (' illegal' if lab['illegal'] else ''), 'INV same' if i0 == i1 else 'INV CHANGES ' + short(i0) + ' -> ' + short(i1)))
        x = c
    return x

def replay(rec, root, ppath, step, cpath, rroot, rpath, pk0, out):
    """one recorded join; returns summary tuple."""
    rows = []; r = nk.LinearNakayamaAlgebra.fromClassName(root)
    P = walk(r, ppath, pk0, 'parent path', rows)
    X = walk(P, step, pk0, 'STEP', rows) if P is not None else None
    MX = walk(X, cpath, pk0, 'child path', rows) if X is not None else None
    R2 = nk.LinearNakayamaAlgebra.fromClassName(rroot) if rroot != 'PARENT' else P
    MR = walk(R2, rpath, pk0, 'root path', rows) if R2 is not None else None
    meet = MX is not None and MR is not None and iso(MX, MR)
    nJ = sum(1 for t in rows if t[2] == 'J!=0'); nNT = sum(1 for t in rows if t[3] == 'NOT-tp')
    ok = all(t[2] != 'GATE REFUSES' for t in rows); inv_ch = [t for t in rows if len(t) > 5 and t[5] != 'INV same']
    ends = (INV(r) == INV(R2)) if R2 is not None else None
    print('\n== %s  root %s  parent path %s  step %s  child path %s || %s path %s' % (rec, root, ppath, step, cpath, rroot, rpath))
    for t in rows: print('   %-12s %-4s %s' % (t[0], t[1], '  '.join(t[2:])))
    print('   replay ok %s | meet iso %s | J != 0 steps %d (%s) | tiltingPlus-false %d | INV changes %d | INV(root) == INV(root\') %s | INV %s'
          % (ok, meet, nJ, ', '.join('%s %s' % (t[0], t[1]) for t in rows if t[2] == 'J!=0'), nNT, len(inv_ch), ends, short(INV(r))), flush=True)
    out.append((rec, ok, meet, nJ, len(inv_ch), len(rows)))
    return P, X

if MODE == 'records':
    t0 = time.time(); pk0 = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
    assert pk0 == search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3304000')) == search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('2400230'))
    A = 'F1 F7 F6 F5 F4'
    RECS = [  # (record, root, parent path, step, child path, root of the other side, its path)  -- sources in the submission
        ('E-233', '3060000', 'F1 R9 R8 R8', 'R8', 'R1 R3', '3060000', A),
        ('E-233 X~P', '3060000', 'F1 R9 R8 R8', 'R8', 'F8 F8', 'PARENT', 'F8 R7 R6 R5'),
        ('E-234 1932', '3060000', 'F1 R8 R9 R9 R8', 'R8', 'F9 R1 R3', '3060000', A),
        ('E-234 1970', '3060000', 'F1 R9 F1 R8 R8', 'R8', 'R1 R1 R3', '3060000', A),
        ('E-234 2003', '3060000', 'F1 R9 F3 R8 R8', 'R8', 'R3 R1 R3', '3060000', A),
        ('E-234 4551', '3060000', 'R8 R9 R9 F1 R8', 'R8', 'R1 F9 R3', '3060000', A),
        ('E-234 4722', '3060000', 'R9 R8 R8 F1 F8', 'R1', 'R3 R8', '3060000', A),
        ('E-239 2026', '3304000', 'F2 R8 R9 R7 R7', 'R7', 'R2', '3304000', 'F2 F1 F6 F5 F2'),
        ('E-244 4606', '3304000', 'F2 F1 R8 R9 R7 R7', 'R7', 'R1 R2 R4', '3304000', 'F2 F1 F3 F6 F5'),
        ('E-244 8063', '3304000', 'F2 R8 R7 R9 F7 R9', 'R9', 'R2 R4 R5', '3304000', 'F2 F1 F6'),
        ('E-244 8170', '3304000', 'F2 R8 R9 F2 R7 R7', 'R7', 'R2 R2 R4', '3304000', 'F2 F1 F6 F5'),
        ('E-244 8196', '3304000', 'F2 R8 R9 F4 R7 R7', 'R7', 'R4 R2 R4', '3304000', 'F2 F1 F6 F5'),
        ('E-244 11123', '3304000', 'R8 R9 R7 R7 F2 F7', 'R2', 'R4 R7 R5', '3304000', 'F2 F1 F6'),
        ('E-254', '2400230', 'F1 F1 F2 F3 R9 R7 R9 F3', 'F3', 'F3 F7 F9', '2400230', 'F1 R4 R5 R9 R7 R9')]
    out = []; kids = {}
    for rec in RECS:
        P, X = replay(*rec, pk0, out); kids[rec[0]] = X
    # E-255: 8427's child is isomorphic to 8063's child (no own join path recorded)
    rows = []; P = walk(nk.LinearNakayamaAlgebra.fromClassName('3304000'), 'F2 R9 R8 R7 R9 R7', pk0, 'parent path', rows); X = walk(P, 'R7', pk0, 'STEP', rows)
    print('\n== E-255 8427: parent path F2 R9 R8 R7 R9 R7 step R7'); [print('   %-12s %-4s %s' % (t[0], t[1], '  '.join(t[2:]))) for t in rows]
    print('   child iso to E-244 8063 child: %s; INV equal %s' % (iso(X, kids['E-244 8063']), INV(X) == INV(kids['E-244 8063'])))
    ia, ib = INV(nk.LinearNakayamaAlgebra.fromClassName('3060000')), INV(nk.LinearNakayamaAlgebra.fromClassName('3304000'))
    print('\nINV separates the F-010 halves: %s (SNF equal %s, signature equal %s, q_S histogram %s vs %s)' % (ia != ib, ia[0] == ib[0], ia[1] == ib[1], ia[2], ib[2]))
    # planted wrong merges: E-233's child path against 3304000's E-239 root path, E-239's child path against 3060000's root path (must fail)
    pl = []
    replay('PLANTED a', '3060000', 'F1 R9 R8 R8', 'R8', 'R1 R3', '3304000', 'F2 F1 F6 F5 F2', pk0, pl)
    replay('PLANTED b', '3304000', 'F2 R8 R9 R7 R7', 'R7', 'R2', '3060000', A, pk0, pl)
    print('planted cross merges: meet iso %s (both must be False); one INV comparison (both plants pair the same two roots): INV(3060000) == INV(3304000) %s; INV changes at their steps %s'
          % ([o[2] for o in pl], INV(nk.LinearNakayamaAlgebra.fromClassName('3060000')) == INV(nk.LinearNakayamaAlgebra.fromClassName('3304000')), [o[4] for o in pl]))
    print('\nSUMMARY (record, replay ok, meet iso, J != 0 steps, INV changes, steps):')
    for o in out: print('  ', o)
    print('merges with >= 1 J != 0 step: %d of %d; with exactly 1: %d; with an INV change: %d; time %.0fs'
          % (sum(o[3] > 0 for o in out), len(out), sum(o[3] == 1 for o in out), sum(o[4] > 0 for o in out), time.time() - t0))

elif MODE == 'mkpkl':
    # the five E-234 joins whose paths were read only from the terminal (1275 1433 1746 1818 1915): write parent pickles in the format
    # skeptic_join.py `pkl:PATH:child` reads, so their joins can be regenerated (paths from skeptic_jcount.py STEPS, E-316)
    import pickle, os; outdir = ARGV[2]; os.makedirs(outdir, exist_ok=True)
    pk0 = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
    for idx, path, step in (('1275', 'F1 F4 R9 R8 R8', 'R8'), ('1433', 'F1 F5 R9 R8 R8', 'R8'), ('1746', 'F1 R5 R9 R8 R8', 'R8'),
                            ('1818', 'F1 R6 R9 R8 R8', 'R8'), ('1915', 'F1 R8 R9 F8 R9', 'R9')):
        rows = []; P = walk(nk.LinearNakayamaAlgebra.fromClassName('3060000'), path, pk0, 'parent path', rows)
        (kind, v), = PP(step); pickle.dump([dict(parent=P, kind=kind, v=v)], open(os.path.join(outdir, 'k%s.pkl' % idx), 'wb'))
        print(idx, 'parent path J', [t[2] for t in rows], 'written')
elif MODE == 'extra':
    # regenerated round 123 by skeptic_join.py side R 3060000 5 / side C pkl:k<i>.pkl:child 3 / meet (totals 8 8 7 7 7 = E-234's table)
    pk0 = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000')); out = []; t0 = time.time()
    EXTRA = [('E-234 1275', 'F1 F4 R9 R8 R8', 'R8', 'F8 F8 F8', 'F1 F5 R9 R8 R7'), ('E-234 1433', 'F1 F5 R9 R8 R8', 'R8', 'F8 F8 F8', 'F1 F6 R9 R8 R7'),
             ('E-234 1746', 'F1 R5 R9 R8 R8', 'R8', 'R1 R3', 'F1 F4 F7 F6 F5'), ('E-234 1818', 'F1 R6 R9 R8 R8', 'R8', 'R1 R3', 'F1 F5 F7 F6 F4'),
             ('E-234 1915', 'F1 R8 R9 F8 R9', 'R9', 'R1 R3', 'F1 F7 F6 F5 F4')]
    for rec, pp, st, cp, rp in EXTRA: replay(rec, '3060000', pp, st, cp, '3060000', rp, pk0, out)
    print('\nSUMMARY', out); print('merges with >= 1 J != 0 step: %d of %d; exactly 1: %d; INV change: %d; time %.0fs'
          % (sum(o[3] > 0 for o in out), len(out), sum(o[3] == 1 for o in out), sum(o[4] > 0 for o in out), time.time() - t0))
elif MODE == 'f037':
    # F-037 / E-032: the library's F-only walk joins 34504030 to 50505000 by 7 right mutations (5 recorded paths, FINDINGS F-037)
    r = nk.LinearNakayamaAlgebra.fromClassName('34504030'); pk0 = search._coxeterKeyOrNone(r); L = nk.LinearNakayamaAlgebra.fromClassName('50505000'); t0 = time.time()
    for p in ('4 2 1 1 2 2 4', '4 1 2 2 2 4 1', '4 1 2 2 2 1 4', '4 1 2 2 1 2 4', '4 1 1 2 2 2 4'):
        rows = []; e = walk(r, ' '.join('F' + v for v in p.split()), pk0, 'F-037', rows)
        print('F-037 path %s: %s | end iso 50505000 %s, iso its dual %s | INV(ends) equal %s' % (p, ' '.join('%s/%s/%s' % (t[2], t[3], 'INV=' if t[5] == 'INV same' else 'INV!') for t in rows),
              iso(e, L), iso(e, pathAlgebra.dualPathAlgebra(L)), INV(r) == INV(L)), flush=True)
    # F-054 (FINDINGS l. 60-90): 30330300 by [7] meets 30330400 by [8, 7]; mixed certificate [7, -9, -10] (negative = left = R here);
    # and the C(2,4,5) pair: 05040330 by [4, 3, 7, 2] meets 60504030 by [-10, -7, -2] (a negative path runs on the opposite algebra = R moves)
    out = []
    for rec, a, pa, b, pb in (('F-054 meet', '30330300', 'F7', '30330400', 'F8 F7'), ('F-054 mixed', '30330300', 'F7 R9 R10', '30330400', ''),
                              ('C(2,4,5) meet', '05040330', 'F4 F3 F7 F2', '60504030', 'R10 R7 R2')):
        pk0 = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName(a)); replay(rec, a, '', '', pa, b, pb, pk0, out)
    print('SUMMARY', out, 'time %.0fs' % (time.time() - t0))
elif MODE == 'control':
    t0 = time.time(); r = nk.LinearNakayamaAlgebra.fromClassName('003033030'); pk0 = search._coxeterKeyOrNone(r); out = []
    # E-318: 003033030 -F3 R8 F3-> P -F3-> X -F6-> 030330030 ; the meet is the LNA itself (zero-length path from 030330030)
    rows = []; MX = walk(r, 'F3 R8 F3 F3 F6', pk0, 'E-318 walk', rows); L = nk.LinearNakayamaAlgebra.fromClassName('030330030')
    print('== control E-318 5-step walk from 003033030')
    for t in rows: print('   %-12s %-4s %s' % (t[0], t[1], '  '.join(t[2:])))
    print('   end iso 030330030: %s | J != 0 steps %d | INV changes %d | INV(003033030) == INV(030330030): %s  FLAGGED: %s  time %.0fs'
          % (iso(MX, L), sum(t[2] == 'J!=0' for t in rows), sum(t[5] != 'INV same' for t in rows), INV(r) == INV(L),
             any(t[2] == 'J!=0' for t in rows) and any(t[5] != 'INV same' for t in rows), time.time() - t0))
