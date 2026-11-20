"""Round 133 (toolsmith), T10 agenda item 1: one counted pass of `refuseNonTilting` over walks that can hold d >= 2 steps.

At every gate-admitted step: d = max_{i != v} dim e_i A e_v (paths i ~> v modulo I), `isTiltingAtVertex`, and -- for steps with
d >= 2 or not tilting -- legality and whether the child keeps the root's Coxeter key.  Non-vacuous step = admitted, legal, key-kept,
d >= 2 (E-338: the only steps where the option can refuse and change the walk).  Pre-registration: workshop/rounds/133/toolsmith.md.

modes (repository root):
  control                      E-318's planted step: tilting must be False, d = 2; option-on depth-1 walk from P must not reach the child
  lib N D ROOT [dual|rdual]    library walk (defaults, fingerprint.Visited) counted, then the option-on walk; LNAs lost/gained with SNF class
                               (dual = pathAlgebra.dualPathAlgebra, as maverick_libwalk.py; rdual = the LNA's relation dual, as merges.py)
  members N                    print merges.py's leftover members at N (one per line)
  nodes N D ROOT [dual|rdual]  post hoc: all nodes (canonicalKey) default vs option on, SNF class of lost nodes
  paths                        recorded mixed F/R paths: E-316 joins (E-321 replay lists), E-255, E-318 walk, E-332 10-step walk
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py lib 11 6 303500030
"""
import sys, os, json, time, copy, collections
sys.path.insert(0, '.')
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import search, mutation, procedure, reduction, pathAlgebra, fingerprint, invariants, lines
from quivermutation import nakayama as nk, arrowPaths as ap, lnaMoves as lm

def maxDim(alg, v):
    q = alg.quiver; rel = procedure.relationsFrom(alg); d = 0
    for i in q.nodes:
        if i == v: continue
        c = ap.allPathsBetween(q, i, v)
        if c: d = max(d, len(c) - len(ap.idealBasis(q, rel, i, v)))
    return d

def inv(A):
    C = invariants.integerCartanMatrix(A); n = len(C)
    S = [[int(C[i][j]) + int(C[j][i]) for j in range(n)] for i in range(n)]
    M = smith_normal_form(sp.Matrix(S), domain=sp.ZZ); d = tuple(sorted(abs(int(M[i, i])) for i in range(n)))
    ev = np.linalg.eigvalsh(np.array(S, dtype=float)); return d, (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))

def childOf(a, v):
    """(legal, child) of right mutation of `a` at v, as the library walk builds it."""
    raw = mutation.quiverMutationAtVertex(copy.deepcopy(a), v)
    if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)): return False, None
    return True, reduction.reducePathAlgebra(raw)

_gate, _tilt = mutation.mutationIsPossibleAtVertex, mutation.isTiltingAtVertex
ST = dict(on=False, base=None, C=collections.Counter(), rows=[])
def countedGate(alg, v):
    ok = _gate(alg, v)
    if ok and ST['on']:
        C = ST['C']; C['admitted'] += 1
        d = maxDim(alg, v); t = _tilt(alg, v); C['d>=2' if d >= 2 else 'd<=1'] += 1
        if not t: C['notTilting'] += 1
        if d >= 2 or not t:
            legal, c = childOf(alg, v)
            kept = legal and search._coxeterKeyOrNone(c) == ST['base']
            C['d%d tilt=%s legal=%s kept=%s' % (min(d, 9), t, legal, kept)] += 1
            if kept and d >= 2: C['NONVACUOUS'] += 1
            if not t and d <= 1: C['ALARM notTilting with d<=1 (contradicts E-338)'] += 1
            if not t and kept: ST['rows'].append(dict(d=d, v=v))
    return ok
mutation.mutationIsPossibleAtVertex = countedGate

def libwalk(start, N, D, refuse, count):
    ST['on'] = count; ST['C'] = collections.Counter(); ST['rows'] = []; ST['base'] = search._coxeterKeyOrNone(start)
    reached = {}
    def visit(A, path):
        if len(A.quiver.edges) != N - 1 or len(A.quiver.nodes) != N: return
        row = lm.asRelLengths(lm._copy(A), N)
        if row is not None:
            name = lines.className(tuple(row))
            if name not in reached or len(path) < len(reached[name]): reached[name] = list(path)
    t0 = time.time()
    search.mutationSearchDepthFirst(copy.deepcopy(start), D, printOutput=False, visitor=visit, visited=fingerprint.Visited(), refuseNonTilting=refuse)
    ST['on'] = False
    return reached, dict(ST['C']), list(ST['rows']), round(time.time() - t0, 1)

def startOf(N, root, side):
    r = nk.LinearNakayamaAlgebra.fromClassName(root)
    if side == 'dual': return pathAlgebra.dualPathAlgebra(r)
    if side == 'rdual': return r.relationDual()
    return r

def lib(N, D, root, side):
    s = startOf(N, root, side); i0 = inv(s)
    off, C, rows, t_off = libwalk(s, N, D, False, True)
    on, _, _, t_on = libwalk(s, N, D, True, False)
    lost, gained = sorted(set(off) - set(on)), sorted(set(on) - set(off))
    lab = lambda nm: 'in' if inv(nk.LinearNakayamaAlgebra.fromClassName(nm)) == i0 else 'CROSS'
    rec = dict(N=N, D=D, root=root, side=side, snfTail=i0[0][-3:], lnasOff=len(off), lnasOn=len(on),
               lost=[(nm, lab(nm), ' '.join(map(str, off[nm]))) for nm in lost], gained=gained,
               crossReachedOff=[nm for nm in off if lab(nm) == 'CROSS'], counts=C, refusableKeptSteps=rows[:20],
               secsOff=t_off, secsOn=t_on)
    print(json.dumps(rec), flush=True)

def nodes(N, D, root, side):
    """post hoc (not pre-registered): every node (canonicalKey), default vs option on; SNF class of the nodes the option loses."""
    s = startOf(N, root, side); i0 = inv(s); seen = {}
    for refuse in (False, True):
        got = {}
        def visit(A, path): got.setdefault(fingerprint.canonicalKey(A), (copy.deepcopy(A), list(path)))
        search.mutationSearchDepthFirst(copy.deepcopy(s), D, printOutput=False, visitor=visit, visited=fingerprint.Visited(), refuseNonTilting=refuse)
        seen[refuse] = got
    lost = [k for k in seen[False] if k not in seen[True]]; gained = [k for k in seen[True] if k not in seen[False]]
    lab = collections.Counter('in' if inv(seen[False][k][0]) == i0 else 'CROSS' for k in lost)
    cross_off = sum(inv(A) != i0 for A, _ in seen[False].values())
    print(json.dumps(dict(N=N, D=D, root=root, side=side, nodesOff=len(seen[False]), nodesOn=len(seen[True]), lost=len(lost),
                          lostByClass=dict(lab), gained=len(gained), crossNodesReachedOff=cross_off,
                          lostPaths=[' '.join(map(str, seen[False][k][1])) for k in lost][:10])), flush=True)

def control():
    sys.path.insert(0, 'tests'); import test_search_soundness as T
    P = T._e318Parent(); base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('003033030'))
    legal, c = childOf(P, 3); ck = fingerprint.canonicalKey(c)
    print('E-318 parent, vertex 3: gate %s, tilting %s, d = %d, legal %s, key kept %s'
          % (_gate(P, 3), _tilt(P, 3), maxDim(P, 3), legal, search._coxeterKeyOrNone(c) == base))
    for refuse in (False, True):
        keys = []
        search.mutationSearchDepthFirst(copy.deepcopy(P), 1, printOutput=False, baseKey=base, refuseNonTilting=refuse,
                                        visitor=lambda A, p: keys.append(fingerprint.canonicalKey(A)))
        print('  depth-1 library walk from P, refuseNonTilting=%s: %d nodes, reaches the E-318 child: %s' % (refuse, len(keys), ck in keys))

def members(N):
    import merges
    class P:   # serial stand-in for the pool
        map = staticmethod(lambda f, xs: list(map(f, xs)))
    _, groups = merges.leftoverGroups(N, P)
    for poly, orbits in groups.items():
        for o in orbits:
            for m in o['members']: print(lines.className(tuple(m)), o['name'])

def walkPath(x, path, base, tag, out):
    for t in path.split():
        kind, v = t[0], int(t[1:]); a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        if not _gate(a, v): out.append((tag, t, 'GATE REFUSES')); return None
        d, tl = maxDim(a, v), _tilt(a, v); legal, c = childOf(a, v)
        if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
        kept = search._coxeterKeyOrNone(c) == base
        out.append((tag, t, d, tl, legal, kept)); x = c
    return x

def paths():
    A = 'F1 F7 F6 F5 F4'
    RECS = [('E-233', '3060000', 'F1 R9 R8 R8', 'R8', 'R1 R3', '3060000', A),
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
        ('E-254', '2400230', 'F1 F1 F2 F3 R9 R7 R9 F3', 'F3', 'F3 F7 F9', '2400230', 'F1 R4 R5 R9 R7 R9'),
        ('E-234 1275', '3060000', 'F1 F4 R9 R8 R8', 'R8', 'F8 F8 F8', '3060000', 'F1 F5 R9 R8 R7'),
        ('E-234 1433', '3060000', 'F1 F5 R9 R8 R8', 'R8', 'F8 F8 F8', '3060000', 'F1 F6 R9 R8 R7'),
        ('E-234 1746', '3060000', 'F1 R5 R9 R8 R8', 'R8', 'R1 R3', '3060000', 'F1 F4 F7 F6 F5'),
        ('E-234 1818', '3060000', 'F1 R6 R9 R8 R8', 'R8', 'R1 R3', '3060000', 'F1 F5 F7 F6 F4'),
        ('E-234 1915', '3060000', 'F1 R8 R9 F8 R9', 'R9', 'R1 R3', '3060000', 'F1 F7 F6 F5 F4'),
        ('E-255 8427', '3304000', 'F2 R9 R8 R7 R9 R7', 'R7', '', None, '')]
    tot = collections.Counter(); t0 = time.time()
    for rec, root, pp, st, cp, rroot, rp in RECS:
        r = nk.LinearNakayamaAlgebra.fromClassName(root); base = search._coxeterKeyOrNone(r); out = []
        P = walkPath(r, pp, base, 'parent', out); X = walkPath(P, st, base, 'STEP', out) if P is not None else None
        if X is not None and cp: walkPath(X, cp, base, 'child', out)
        if rroot: walkPath(P if rroot == 'PARENT' else nk.LinearNakayamaAlgebra.fromClassName(rroot), rp, base, 'rootpath', out)
        ends = 'same LNA (%s)' % root if rroot in (root, 'PARENT') else 'E-255: no own join path' if rroot is None else 'DISTINCT'
        tally(out, tot)
        print('%-12s ends %-22s steps %2d | d>=2 %d | non-vacuous %d | refused %d (%s) | refused d<=1 %d'
              % (rec, ends, len(out), sum(o[2] >= 2 for o in out if len(o) > 3), sum(o[2] >= 2 and o[4] and o[5] for o in out if len(o) > 3),
                 sum(not o[3] for o in out if len(o) > 3), ','.join('%s %s d%d' % (o[0], o[1], o[2]) for o in out if len(o) > 3 and not o[3]),
                 sum(not o[3] and o[2] <= 1 for o in out if len(o) > 3)), flush=True)
    for label, root, side, path in (('E-318 walk (control)', '003033030', 'root', 'F3 R8 F3 F3 F6'),
                                    ('E-332 walk', '303305000', 'dual', 'F6 F11 F7 F11 F7 F7 F3 F3 F4 F7')):
        s = startOf(0, root, side); base = search._coxeterKeyOrNone(s); out = []; walkPath(s, path, base, 'walk', out); tally(out, tot)
        print('%-20s %s' % (label, ' '.join('%s:d%d%s%s' % (o[1], o[2], '' if o[3] else '/REFUSED', '' if o[5] else '/moved') if len(o) > 3 else o[1] + ':' + o[2] for o in out)), flush=True)
    print('TOTAL', dict(tot), '%.0fs' % (time.time() - t0))

def tally(out, tot):
    for o in out:
        if len(o) <= 3: tot['gate refuses'] += 1; continue
        _, _, d, t, legal, kept = o; tot['steps'] += 1
        tot['d>=2' if d >= 2 else 'd<=1'] += 1; tot['d>=2 nonvacuous' if d >= 2 and legal and kept else 'x'] += 0
        if d >= 2 and legal and kept: tot['nonvacuous'] += 1; tot['nonvacuous tilting=%s' % t] += 1
        if not t: tot['refused d=%d kept=%s' % (d, kept)] += 1

if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == 'control': control()
    elif a[0] == 'lib': lib(int(a[1]), int(a[2]), a[3], a[4] if len(a) > 4 else 'root')
    elif a[0] == 'members': members(int(a[1]))
    elif a[0] == 'nodes': nodes(int(a[1]), int(a[2]), a[3], a[4] if len(a) > 4 else 'root')
    elif a[0] == 'paths': paths()
