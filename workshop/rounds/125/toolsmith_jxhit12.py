"""Round 125 (toolsmith): independent checks of the n = 12 crossing found by toolsmith_jxclass.py from 0040033030 (key group of 47, E-314 /
E-317: 14 LNAs SNF tail (1,1,20), 33 LNAs (1,1,12)).  Rebuild P = 0040033030 . F3 R9 F3, take F3 (and the dual-side mirror: dual . F9 R3 R3,
then R3), then with round-122 code that does not use the library's Cartan routine: brute-force Cartan == integerCartanMatrix, det C, acyclic,
SNF(S) by sympy and the determinantal divisor d_11 of S (rank 11, so d_11 = product of the nonzero invariant factors), gate / perI /
isTiltingAtVertex / key at the step; then a J = 0 key-kept F/R BFS from the child, radius R (default 2), listing LNAs reached with their SNF tail.
usage (repo root): R=2 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxhit12.py"""
import sys, os, time, copy, collections, itertools, math
sys.path.insert(0, '.'); sys.argv = ['x', 'none']
_s = open('workshop/rounds/122/toolsmith_jxsearch.py').read(); exec(compile(_s.split("if __name__ == '__main__':")[0], 'jx', 'exec'))
_h = open('workshop/rounds/122/toolsmith_jxhit.py').read(); exec(compile(_h[_h.index('def step('):_h.index('lab = labels(N)')], 'jxhit', 'exec'))
N = 12; ROOT = '0040033030'; R = int(os.environ.get('R', '2'))
r = nk.LinearNakayamaAlgebra.fromClassName(ROOT); pk0 = search._coxeterKeyOrNone(r); t0 = time.time()
snft = lambda C: INV(C)[0][-3:]
print('root', ROOT, 'SNF tail', snft(ct.lnaCartanMatrix(N, tuple(int(c) for c in ROOT))), '; 33-class root 0400033030 tail', snft(ct.lnaCartanMatrix(N, (0, 4, 0, 0, 0, 3, 3, 0, 3, 0))))
kids = []
for side, path, st in (('root', 'F3 R9 F3', 'F3'), ('dual', 'F9 R3 R3', 'R3')):
    x = copy.deepcopy(r) if side == 'root' else pathAlgebra.dualPathAlgebra(r); js = []
    for mv in path.split():
        x, J = step(x, mv); js.append(J); assert search._coxeterKeyOrNone(x) == pk0
    kind, v = st[0], int(st[1:]); a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    c, J = step(x, st); kids.append(c)
    print('\n%s %s then %s: path J != 0? %s; step: gate yes, J != 0 %s, isTiltingAtVertex %s, key kept %s' % (side, path, st, js, J, mutation.isTiltingAtVertex(a, v), search._coxeterKeyOrNone(c) == pk0))
    for nm, A in (('parent', x), ('child', c)):
        C = invariants.integerCartanMatrix(A); B = brute_cartan(A); S = [[C[i][j] + C[j][i] for j in range(N)] for i in range(N)]
        print('  %s: acyclic %s, brute-force Cartan == library %s, det C %s, SNF tail %s, d_11(S) = %s, rank S %d' % (nm, nx.is_directed_acyclic_graph(A.quiver), B == C,
              sp.Matrix(C).det(), snft(C), detdiv(S, N - 1), sp.Matrix(S).rank()), flush=True)
print('children isomorphic:', iso(kids[0], kids[1]), ' child0 ~ dual(child1):', iso(kids[0], pathAlgebra.dualPathAlgebra(kids[1])))
c = kids[0]; seen = {fingerprint.canonicalKey(c)}; fr = [c]; found = collections.Counter()
for d in range(1, R + 1):
    nf = []
    for x in fr:
        if list(nx.simple_cycles(x.quiver)): continue
        for kind, v in [(kd, v) for kd in 'FR' for v in sorted(x.quiver.nodes)]:
            a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
            if not mutation.mutationIsPossibleAtVertex(a, v) or perI(a, v): continue
            raw = mutation.quiverMutationAtVertex(a, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if kind == 'R': ch = pathAlgebra.dualPathAlgebra(ch)
            if search._coxeterKeyOrNone(ch) != pk0: continue
            k = fingerprint.canonicalKey(ch)
            if k is None or k in seen: continue
            seen.add(k); nf.append(ch); nm = lna_name(ch, N)
            if nm: found[(nm, snft(invariants.integerCartanMatrix(ch)), d, kind + str(v))] += 1
    fr = nf
    print('child J = 0 ball radius %d: %d nodes, LNAs (name, SNF tail, radius, last step) %s  %.0fs' % (d, len(seen), sorted(found), time.time() - t0), flush=True)
