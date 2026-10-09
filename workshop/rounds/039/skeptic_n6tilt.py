"""Round 039 (skeptic): replay an LNA -> M and a hit -> M path of E-136 (n = 6, class 0) step by step.
Each step: gate (mutationIsPossibleAtVertex), tiltingPlus (Ladkani 2.3(c)), Cartan congruence (Ladkani 3.6), Coxeter key.
Usage: skeptic_n6path.py LNA_SECONDS HIT_SECONDS [HITS=0,1,2]"""
import sys, time
lnaS = float(sys.argv[1]); hitS = float(sys.argv[2]); which = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '0').split(',')]
src = open('workshop/rounds/038/toolsmith_n6close.py').read().split("budget = args")[0]
sys.argv = ['x', '--plan', '1']; exec(compile(src, 'c', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]
h15 = h15.replace("from quivermutation import nakayama as nk", "from quivermutation import nakayama as nk")
exec(compile(h15, 'h015', 'exec'))
key0 = order[0]
verts = list(range(1, 7))
TILT = False
def bfs(starts, secs):
    seen = {}; fr = []; t0 = time.time(); nonek = 0
    def add(x, par, v, depth):
        nonlocal nonek
        ck = fingerprint.canonicalKey(x)
        if ck is None: nonek += 1; return
        if ck in seen: return
        seen[ck] = (par, v, x, depth); fr.append((x, ck, depth))
    for x in starts: add(x, None, None, 0)
    done = False
    while fr and not done:
        cur = fr[:]; fr.clear()
        for j, (alg, ck0, dp) in enumerate(cur):
            if time.time() - t0 > secs: done = True; break
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                if TILT and not tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v): continue
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == key0: add(ch, ck0, v, dp + 1)
    return seen, nonek
def path(seen, ck):
    out = []
    while seen[ck][0] is not None:
        par, v, x, d = seen[ck]; out.append((par, v, ck)); ck = par
    return ck, out[::-1]
def replay(seen, ck, label):
    root, steps = path(seen, ck)
    print('--- %s: %d steps from root key %s' % (label, len(steps), str(root)[:60]))
    bad = 0
    for par, v, c in steps:
        A = seen[par][2]; B = seen[c][2]; rels = procedure.relationsFrom(A)
        gate = mutation.mutationIsPossibleAtVertex(A, v); tp = tiltingPlus(A.quiver, rels, v)
        raw = mutation.quiverMutationAtVertex(A, v); ch = reduction.reducePathAlgebra(raw)
        same = fingerprint.canonicalKey(ch) == c
        cong = len(ch.vertices()) == 6 and bool((rplus(A, v, verts).dot(cartan(A)).dot(rplus(A, v, verts).T) == cartan(ch)).all())
        kd = kerdims(A, v, rels); kd = {i: x for i, x in kd.items() if x}
        ok = gate and tp and cong and same
        bad += (not ok)
        print('  v=%d gate=%s tiltingPlus=%s cartanCong=%s recompute-same=%s coxkey-pres=%s ker=%s |A|dim=%d' % (v, gate, tp, cong, same, search._coxeterKeyOrNone(ch) == key0, kd, int(sum(cartan(A).flatten()))))
    print('  steps failing any test:', bad); return bad
def kerdims(alg, k, rels):
    quiver = alg.quiver; outs = ap.arrowsOutOf(quiver, k); res = {}
    for i in quiver.nodes:
        if i == k: continue
        P = ap.allPathsBetween(quiver, i, k)
        if not P: res[i] = 0; continue
        dimV = len(P) - len(ap.idealBasis(quiver, rels, i, k)); rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, v in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): row[(b, kk)] = v
            rows.append(row)
        res[i] = dimV - rank(rows)
    return res
TILT = True
S1, n1 = bfs(classes[key0], lnaS); print('TILTING-ONLY LNA side', len(S1), flush=True)
starts = [build([tuple(x) for x in h[0]], h[1], list(range(1, 7))) for h in hits]
tot = 0
for j in range(len(starts)):
    S2, n2 = bfs([starts[j]], hitS); sh = set(S1) & set(S2)
    print('hit %d (v=%d): tilting-only seen %d shared with tilting-only LNA side %d' % (j, hits[j][3], len(S2), len(sh)), flush=True)
    if sh and j in which:
        m = min(sh, key=lambda k: S1[k][3] + S2[k][3]); tot += replay(S1, m, 'LNA -> M'); tot += replay(S2, m, 'hit -> M')
print('TOTAL failing steps', tot)
