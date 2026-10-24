"""Round 085 (skeptic) response, item 3: Hom-verify the lnaMoves edges.  Every move edge u -> v (one rule application, sequence of signed vertices)
is replayed step by step with the round-079 `replay` check (Hom(T,T[-1]) = Hom(T,T[1]) = 0 and End(T) dims == Cartan(child)^T for each mutation, S = -1),
and the final algebra must standardise to the predicted class name.
usage: skeptic_hom.py comp            all edges of the 7-row component of 3304000 (n = 9)
       skeptic_hom.py all                    every move edge of every row (n = 9), 480 s budget
       skeptic_hom.py sample SEED K   K random edges (one per row) (n = 9)
Run from the repository root."""
import sys, random, time
A_ = sys.argv; MODE = A_[1]
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import numpy as np
from quivermutation import invariants, lnaMoves, lines
n = 9
def replaySeq(x0, seq):
    ok = True; x = x0
    for s in seq:
        kind, v = ('F', s) if s > 0 else ('R', -s)
        a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        Vs, res = stepTest(a, v, -1)
        zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
        c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
        C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
        ok = ok and zero and bool((np.array(res[0], dtype=int) == C.T).all())
        x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return ok, x
def edges(row):
    s = ''.join(map(str, row)); out = []
    for name, seq in lnaMoves.movesByRule(n, [int(c) for c in s]).items(): out.append((s, name, list(seq)))
    return out
def check(s, name, seq):
    x0 = nk.LinearNakayamaAlgebra(n, [int(c) for c in s])
    ok, x = replaySeq(x0, seq)
    _, rl, _ = lnaMoves.standardise(x)
    got = lines.className(rl) if rl is not None else None
    return ok, got == name, got
if MODE == 'comp':
    comp = ['3304000', '3004000', '3024000', '3004002']
    # full component by undirected closure of the move graph
    seen = set(); stack = ['3304000']; E = []
    allrows = [''.join(map(str, r)) for r in nk.allRelationLengths(n)]
    fw = {}
    for r in allrows:
        for e in edges(r): fw.setdefault(e[0], []).append(e)
    undirected = {}
    for u, es in fw.items():
        for (_, v, _) in es: undirected.setdefault(u, set()).add(v); undirected.setdefault(v, set()).add(u)
    while stack:
        u = stack.pop()
        if u in seen: continue
        seen.add(u); stack += list(undirected.get(u, []))
    print('component of 3304000:', sorted(seen), len(seen))
    tot = bad = 0
    for u in sorted(seen):
        for (s, v, seq) in fw.get(u, []):
            ok, hit, got = check(s, v, seq); tot += 1; bad += not (ok and hit)
            print(s, '->', v, 'seq', seq, 'Hom-verified', ok, 'predicted class reached', hit)
    print('edges', tot, 'failing', bad)
elif MODE == 'all':
    t0 = time.time(); tot = bad = rows = 0
    for r in [''.join(map(str, r)) for r in nk.allRelationLengths(n)]:
        rows += 1
        for (s, v, seq) in edges(r):
            ok, hit, got = check(s, v, seq); tot += 1; bad += not (ok and hit)
            if not (ok and hit): print('FAIL', s, v, seq, ok, hit, got)
        if time.time() - t0 > 480: print('TIME BUDGET HIT'); break
    print('rows', rows, 'edges', tot, 'failing', bad, 'secs %.0f' % (time.time() - t0))
else:
    rnd = random.Random(int(A_[2])); K = int(A_[3])
    allrows = [''.join(map(str, r)) for r in nk.allRelationLengths(n)]
    rnd.shuffle(allrows); pool = []
    t0 = time.time(); tot = bad = 0
    for r in allrows:
        es = edges(r)
        if not es: continue
        s, v, seq = rnd.choice(es)
        ok, hit, got = check(s, v, seq); tot += 1; bad += not (ok and hit)
        if not (ok and hit): print('FAIL', s, v, seq, ok, hit, got)
        if tot >= K or time.time() - t0 > 480: break
    print('sampled edges (one per row, random)', tot, 'failing', bad, 'secs %.0f' % (time.time() - t0))
