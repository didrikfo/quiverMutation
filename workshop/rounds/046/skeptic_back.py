"""Round 046 (skeptic): bounded tilting-only reach test.  From the CHILD of a recorded step (rebuilt from arrows/relations), BFS with gate, legal, vertex-
preserving, J = 0 steps and key kept, up to maxexp expansions; report whether the walk reaches the canonicalKey of an LNA or dual LNA of the class (a hit proves the
child is derived equivalent to the class via tilting steps (Ladkani 2.3(c)); a miss proves nothing).  Controls: kind ctrl (a J = 0 child, in class by construction).
Usage: skeptic_back.py in.pkl kind first last maxexp"""
import sys, pickle, ast, time
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
pk, kind, lo, hi, maxexp = ARGV[1], ARGV[2], int(ARGV[3]), int(ARGV[4]), int(ARGV[5])
def mk(edges, reprr, nodes):
    arrows = [(a, b) for a, b, k in edges]; rels = []
    for d in ast.literal_eval(reprr): rels.append([[p[0][0]] + [a[1] for a in p] for p in d])
    return build(arrows, rels, nodes)
recs = [x for x in pickle.load(open(pk, 'rb')) if x['kind'] == kind and x['depth'] >= (7 if kind == 'fail' else 3)][lo:hi]
n = recs[0]['n']; cls = recs[0]['cls']
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)): classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
targets = {fingerprint.canonicalKey(a) for a in classes[base]}
for idx, r in enumerate(recs, lo):
    try: A = mk(r['pe'], r['pr'], r['V']); ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, r['v'])); perI(ch, r['V'][0])
    except Exception as e: print(idx, 'rebuild failed', repr(e)[:50], flush=True); continue
    seen = {fingerprint.canonicalKey(ch)}; frontier = [ch]; LV = {id(ch): 0}; nexp = 0; hit = None; t0 = time.time()
    if fingerprint.canonicalKey(ch) in targets: hit = 0
    while frontier and nexp < maxexp and hit is None:
        cur, frontier = frontier, []
        for alg in cur:
            if nexp >= maxexp or hit is not None: break
            nexp += 1
            if list(nx.simple_cycles(alg.quiver)): continue
            V = sorted(alg.quiver.nodes); rels = procedure.relationsFrom(alg)
            for v in V:
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                if perI(alg, v): continue
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                c2 = reduction.reducePathAlgebra(raw)
                if sorted(c2.vertices()) != V or search._coxeterKeyOrNone(c2) != base: continue
                k = fingerprint.canonicalKey(c2)
                if k in seen: continue
                seen.add(k); LV[id(c2)] = LV[id(alg)] + 1; frontier.append(c2)
                if k in targets: hit = LV[id(c2)]; break
    print(kind, idx, 'parent depth', r['depth'], 'expanded', nexp, 'closed' if not frontier else 'open', 'HIT at distance %s' % hit if hit is not None else 'no LNA reached', '%.0fs' % (time.time() - t0), flush=True)
