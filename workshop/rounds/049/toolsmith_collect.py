"""Round 049 (toolsmith) copy of rounds/046/skeptic_collect.py that ALSO stores the pickled algebra objects (keys parentObj, childObj), so parallel-arrow parents can be searched from. Original docstring: Round 046 (skeptic), T10 (i): collect key-keeping steps of the E-151 walk (as rounds/045/toolsmith_guardaudit.py mode A:
BFS over canonicalKey-deduplicated algebras, seeds = LNAs and duals of key class cls, expansion only into key-kept children).
Per gate-admitted legal vertex-preserving key-kept step records Cartan matrices of parent and child, J, tiltingPlus, depth, and
(arrows, relations) of both.  Keeps ALL J != 0 key-kept steps (distinct by parent canonicalKey, v) and up to NCTRL J = 0 control steps.
Usage: skeptic_collect.py n cls maxnodes out.pkl [NCTRL]    (run from repo root)"""
import sys, time, pickle, copy
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
n, cls, maxnodes, out = int(_a[1]), int(_a[2]), int(_a[3]), _a[4]; NCTRL = int(_a[5]) if len(_a) > 5 else 200
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
seen = set(); frontier = []; LV = {}
def add(alg, lv=0):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg); LV[id(alg)] = lv
for s in classes[base]: add(s)
recs = []; dist = set(); nexp = 0; t0 = time.time(); nctrl = 0; stepsJ = 0
def ints(alg): return [[int(x) for x in r] for r in cartan(alg).tolist()]
while frontier and nexp < maxnodes:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if nexp >= maxnodes: break
        nexp += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        V = sorted(alg.quiver.nodes); rels = procedure.relationsFrom(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: continue
            if search._coxeterKeyOrNone(ch) != base: continue
            J = perI(alg, v); tp = bool(tiltingPlus(alg.quiver, rels, v))
            ck = fingerprint.canonicalKey(alg)
            if J or not tp:
                stepsJ += 1
                if (ck, v) not in dist:
                    dist.add((ck, v))
                    recs.append(dict(parentObj=copy.deepcopy(alg), childObj=copy.deepcopy(ch), kind='fail', n=n, cls=cls, depth=LV[id(alg)], v=v, V=V, J=J, tp=tp, CA=ints(alg), CB=ints(ch),
                                     pe=sorted(alg.quiver.edges(keys=True)), pr=repr(rels), ce=sorted(ch.quiver.edges(keys=True)), cr=repr(procedure.relationsFrom(ch)), R=[[int(x) for x in r] for r in rplus(alg, v, V).tolist()]))
            elif nctrl < NCTRL and (nexp % 7 == 0):
                nctrl += 1
                recs.append(dict(parentObj=copy.deepcopy(alg), childObj=copy.deepcopy(ch), kind='ctrl', n=n, cls=cls, depth=LV[id(alg)], v=v, V=V, J=J, tp=tp, CA=ints(alg), CB=ints(ch),
                                 pe=sorted(alg.quiver.edges(keys=True)), pr=repr(rels), ce=sorted(ch.quiver.edges(keys=True)), cr=repr(procedure.relationsFrom(ch)), R=[[int(x) for x in r] for r in rplus(alg, v, V).tolist()]))
            add(ch, LV[id(alg)] + 1)
print('n', n, 'cls', cls, 'seeds', len(classes[base]), 'expanded', nexp, 'frontier', len(frontier), 'fail steps', stepsJ, 'distinct fail', len(dist), 'ctrl', nctrl, '%.0fs' % (time.time() - t0))
pickle.dump(recs, open(out, 'wb'))
