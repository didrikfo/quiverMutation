"""Round 033 (experimentalist): hand-built both-die squares (Part 1) and walk reach (Part 2).
Usage: experimentalist_bothdie.py hand | walk n class budget_sec [max_exp]
hand: build squares 1->2->4, 1->3->4 with v=4 out-arrows b1:4->5, b2:4->6; relation p1 b1 = p2 b1 (commutes) and
  monomials killing p1 b2 and p2 b2 (several ways); print gate, kerdim, tiltingPlus, Coxeter key, key in LNA key sets at n=6.
walk: guarded BFS (as rounds/023) of each class at n; for every visited algebra test whether its canonicalKey equals a hand algebra's;
  and tally the 'both-die' pattern (E-122 pattern) per class."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); _a = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = _a
sq = [(1,2),(1,3),(2,4),(3,4),(4,5),(4,6)]
com = [[1,2,4,5],[1,3,4,5]]
H = {
 'M1 monos 2-4-6, 3-4-6':          (sq, [com, [[2,4,6]], [[3,4,6]]]),
 'M2 monos 1-2-4-6, 1-3-4-6':      (sq, [com, [[1,2,4,6]], [[1,3,4,6]]]),
 'M3 mono 1-2-4-6, 3-4-6':         (sq, [com, [[1,2,4,6]], [[3,4,6]]]),
 'M4 comm 1-2-4-6=1-3-4-6 +mono':  (sq, [com, [[1,2,4,6],[1,3,4,6]]]),   # D: commutes both (control, J!=0 gate?)
 'M5 half: mono 3-4-6 only':       (sq, [com, [[3,4,6]]]),
}
def Wcode(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v)
    for R in rels:
        if len(R) < 2: continue
        for b1 in outs:
            if not all(len(q) >= 2 and q[-1] == b1 and q[-2][1] == v for q in R): continue
            b2 = [b for b in outs if b != b1][0]
            if ap.isInIdeal(Q, rels, ap.combination({q[:-1]: c for q, c in R.items()})): continue
            if ap.isInIdeal(Q, rels, ap.combination({q[:-1] + (b2,): c for q, c in R.items()})): return True
    return False
def perI(alg, v):
    """dim J_i = dim ker g_i for each i (i != v)."""
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v); res = {}
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        dimV = len(P) - len(ap.idealBasis(Q, rels, i, v)); rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1])).items(): row[(b, kk)] = x
            rows.append(row)
        d = dimV - (rank(rows) if rows else 0)
        if d: res[i] = d
    return res
def lnaKeys(n):
    ks = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            ks.setdefault(search._coxeterKeyOrNone(alg), 0); ks[search._coxeterKeyOrNone(alg)] += 1
    return ks
mode = _a[1]
if mode == 'hand':
    ks = lnaKeys(6); print('n=6 LNA key classes', len(ks))
    for name, (arr, rels) in H.items():
        A = build(arr, rels, [1,2,3,4,5,6]); key = search._coxeterKeyOrNone(A)
        print(name, '| gate', mutation.mutationIsPossibleAtVertex(A, 4), '| Wcode', Wcode(A, 4), '| J_i', perI(A, 4),
              '| tiltingPlus', bool(tiltingPlus(A.quiver, procedure.relationsFrom(A), 4)), '| key in n=6 LNA keys', key in ks, '| key', key)
        print('    rels', A.rels)
def gen(n, ks):
    """yield (core name, algebra, v, gate, sum dim J_i, coxeter key) for every both-die core + attachments on n vertices."""
    cores = {
     'A(2,2)': ([(1,2),(1,3),(2,4),(3,4),(4,5),(4,6)], [1,2,4], [1,3,4], 4, 5, 6),
     'B(1,2)': ([(1,3),(1,2),(2,3),(3,4),(3,5)],       [1,3],   [1,2,3], 3, 4, 5),
     'C(1,3)': ([(1,4),(1,2),(2,3),(3,4),(4,5),(4,6)], [1,4],   [1,2,3,4], 4, 5, 6),
     'D(2,3)': ([(1,2),(2,4),(1,3),(3,5),(5,4),(4,6),(4,7)], [1,2,4], [1,3,5,4], 4, 6, 7),
    }
    for cname, (arr0, p1, p2, v, b1, b2) in cores.items():
        k0 = len({x for e in arr0 for x in e}); extra = n - k0
        if extra < 0: continue
        com0 = [p1 + [b1], p2 + [b1]]
        def kills(p):
            full = p + [b2]; suf = p[-2:] + [b2]
            return [[full]] if suf == full else [[full], [suf]]
        for ka in kills(p1):
            for kb in kills(p2):
                base_rels = [com0, ka, kb]
                newv = list(range(k0 + 1, k0 + 1 + extra)); nodes0 = list(range(1, k0 + 1))
                def attach(i, nodes, arrows):
                    if i == len(newv): yield arrows; return
                    w = newv[i]
                    for u in nodes:
                        for dirn in (0, 1):
                            yield from attach(i + 1, nodes + [w], arrows + [(u, w) if dirn == 0 else (w, u)])
                for att in attach(0, nodes0, []):
                    arrows = arr0 + att
                    for zmask in range(2 ** len(att) if extra else 1):
                        rels = list(base_rels)
                        for j, (x, y) in enumerate(att):
                            if zmask >> j & 1:
                                for (c, d) in arrows:
                                    if d == x and (c, d) != (x, y): rels.append([[c, x, y]])
                                    if c == y and (c, d) != (x, y): rels.append([[x, y, d]])
                        try:
                            A = build(arrows, rels, list(range(1, n + 1)))
                            if any(ap.isIllegalRelation(A.quiver, rr) for rr in procedure.relationsFrom(A)): yield (cname, None, v, 'illegal', 0, None); continue
                            gate = mutation.mutationIsPossibleAtVertex(A, v)
                            kd = sum(perI(A, v).values()) if gate else 0
                            yield (cname, A, v, gate, kd, search._coxeterKeyOrNone(A))
                        except Exception as e:
                            yield (cname, None, v, 'error ' + type(e).__name__, 0, None)
if mode == 'enum':
    n = int(_a[2]); ks = lnaKeys(n); res = Counter(); ex = {}
    for cname, A, v, gate, kd, key in gen(n, ks):
        if A is None: res[(cname, gate)] += 1; continue
        if not gate: res[(cname, 'gate refuses')] += 1; continue
        res[(cname, 'gate ok', 'J!=0' if kd else 'J=0', 'LNAkey' if key in ks else 'noLNAkey')] += 1
        if kd and key in ks: ex.setdefault(cname, (A.quiver.edges, A.rels))
    print('n', n)
    for k, c in sorted(res.items(), key=str): print(c, k)
    for k, (ar, rl) in ex.items(): print('J!=0 + LNA key example', k, list(ar), rl)
if mode == 'reach':
    # for each both-die (gate ok, J != 0) hand algebra whose key is an LNA key, BFS the matching class(es) and look for its canonicalKey.
    n, budget = int(_a[2]), float(_a[3]); maxexp = int(_a[4]) if len(_a) > 4 else 10**9
    classes = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
        for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
            classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
    targets = {}   # canonicalKey -> (cname, key)
    for cname, A, v, gate, kd, key in gen(n, classes):
        if A is None or not gate or not kd or key not in classes: continue
        ck = fingerprint.canonicalKey(A); targets[ck] = (cname, key)
    print('n', n, 'both-die hand algebras with LNA key (distinct canonicalKeys):', len(targets), Counter(c for c, k in targets.values()))
    bykey = Counter(k for c, k in targets.values())
    t0 = time.time()
    for key, cnt in bykey.items():
        ci = sorted(classes, key=lambda k: (len(classes[k]), str(k))).index(key)
        seen = set(); frontier = []; nalg = 0
        def add(alg):
            k = fingerprint.canonicalKey(alg)
            if k is not None:
                if k in seen: return
                seen.add(k)
            frontier.append(alg)
        for sA in classes[key]: add(sA)
        lvl = 0; closed = False; tb = time.time(); hits = set()
        while frontier and nalg < maxexp and time.time() - tb < budget:
            cur, frontier[:] = frontier[:], []; lvl += 1
            for alg in cur:
                if nalg >= maxexp or time.time() - tb > budget: frontier.extend([alg]); continue
                nalg += 1
                if list(nx.simple_cycles(alg.quiver)): continue
                for v in sorted(alg.vertices()):
                    if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                    raw = mutation.quiverMutationAtVertex(alg, v)
                    if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                    ch = reduction.reducePathAlgebra(raw)
                    if search._coxeterKeyOrNone(ch) == key: add(ch)
        closed = not frontier
        hits = {c for c in targets if targets[c][1] == key and c in seen}
        print('class idx', ci, 'size', len(classes[key]), 'targets', cnt, 'expanded', nalg, 'seen', len(seen), 'levels', lvl, 'CLOSED' if closed else 'capped', 'HITS', len(hits), '%.0fs' % (time.time() - tb), flush=True)
