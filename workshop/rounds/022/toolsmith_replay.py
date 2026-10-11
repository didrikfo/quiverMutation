"""Replay the 10 E-086 key-moved 'M' parents (n = 8 class 2) under the fixed library (round 022).
Usage: toolsmith_replay.py [--old]   (--old: monkeypatch the pre-E-091 reduceAgainstPivots, head-only reduction)
Reads workshop/rounds/014/scholar_walk_n8_c2.txt (lines 'M (depth, rels, v, path)'), rebuilds each parent by
the path from the class's algebra list (as workshop/rounds/015/theorist_cartan.py), then for the step at v:
gate, tiltingPlus, Cartan discrepancy (procedure.cartanDiscrepancy) and key(child) == key(parent)."""
import sys, ast, time
sys.argv, args = ['x'], sys.argv[1:]
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
n, c = 8, 2
if '--old' in args:
    from fractions import Fraction
    def _old(comb, pivots):
        row = {k: Fraction(v) for k, v in comb.items()}
        while row:
            head = min(row)
            if head not in pivots: break
            factor = row[head]; pr = pivots[head]
            row = {k: row.get(k, Fraction(0)) - factor * pr.get(k, Fraction(0)) for k in set(row) | set(pr)}
            row = {k: v for k, v in row.items() if v != 0}
        return row
    ap.reduceAgainstPivots = _old
lines = []
for l in open('workshop/rounds/014/scholar_walk_n8_c2.txt'):
    if l.startswith('M '): lines.append(ast.literal_eval(l[2:].strip()))
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[c]
kept = 0
for d, rels, v, path in lines:
    alg = classes[base][path[0]]
    for w in path[1:]:
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, w))
    rebuilt = [[list(p) for p in r] for r in alg.rels]
    same = (rebuilt == rels)
    pk = search._coxeterKeyOrNone(alg)
    gate = mutation.mutationIsPossibleAtVertex(alg, v)
    tp = tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v)
    rl = procedure.relationsFrom(alg)
    t0 = time.time()
    try:
        nq, nr = procedure.mutateAtVertex(alg.quiver, rl, v, checkCartan=True); err = None
    except procedure.CartanCongruenceError as e:
        err = str(e)[:80]
    t1 = time.time()
    ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(alg, v))
    ck = search._coxeterKeyOrNone(ch)
    keep = ck == pk
    kept += keep
    print('v', v, 'path', path, 'parentRelsMatchFile', same, 'gate', gate, 'tiltingPlus', tp,
          'cartanAssertion', 'PASS' if err is None else 'FAIL ' + err, 'keyKept', keep, 'ms', round(1000 * (t1 - t0)))
print('keyKept in', kept, 'of', len(lines))
