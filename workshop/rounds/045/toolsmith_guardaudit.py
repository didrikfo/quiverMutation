"""Round 045 (toolsmith): T10 guard audit (a).  Cost of an added J = 0 / tiltingPlus check on a key-guarded walk.
Usage: toolsmith_guardaudit.py n cls maxnodes mode          (run from repo root)
  mode 'A'  : walk with the guard as is (key guard only).  For every gate-admitted, legal, vertex-preserving step the three
              costs are timed separately (step = mutate+reduce+key; J = perI on the parent; tp = tiltingPlus on the parent),
              and (J == 0, tiltingPlus, key kept) is tallied.  The walk is the one the library guard gives.
  mode 'J'  : the same walk, but a step is refused (before it is mutated) when J != 0.     [check placed after the gate]
  mode 'T'  : same with tiltingPlus.
  mode 'G'  : as A, but never expands into a child of a tiltingPlus-failing step (tguard, as skeptic_c2).
  A and G also tally depth (BFS level of the parent), distinct (parent,v) and whether the parent was reached through a failing key-keeper.
Walk = BFS over canonicalKey-deduplicated algebras (as rounds/043/skeptic_c2.py), seeds = LNAs and duals of key class `cls`
(classes sorted by (size, str(key))), expansion only into children whose key == class key; stops after `maxnodes` expansions."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
n, cls, maxnodes, mode = int(_a[1]), int(_a[2]), int(_a[3]), _a[4]
pc = time.perf_counter
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
seen = set(); frontier = []; LV = {}; TAINT = {}; FAIL = []   # per-alg BFS level / 'reached through a non-tilting key-keeper'
def add(alg, lv=0, taint=False):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg); LV[id(alg)] = lv; TAINT[id(alg)] = taint
for s in classes[base]: add(s)
T = Counter(); N = Counter(); tab = Counter(); refused = 0; nexp = 0; t0 = pc(); level = 0
while frontier and nexp < maxnodes:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if nexp >= maxnodes: break
        nexp += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        V = sorted(alg.quiver.nodes); rels = procedure.relationsFrom(alg)
        for v in V:
            t = pc(); g = mutation.mutationIsPossibleAtVertex(alg, v); T['gate'] += pc() - t
            if not g: continue
            N['gate'] += 1
            if mode in 'AG':
                t = pc(); J = perI(alg, v); T['J'] += pc() - t; N['J'] += 1
                t = pc(); tp = bool(tiltingPlus(alg.quiver, rels, v)); T['tp'] += pc() - t; N['tp'] += 1
            elif mode == 'J':
                t = pc(); J = perI(alg, v); T['J'] += pc() - t; N['J'] += 1
                if J: refused += 1; continue
            elif mode == 'T':
                t = pc(); tp = bool(tiltingPlus(alg.quiver, rels, v)); T['tp'] += pc() - t; N['tp'] += 1
                if not tp: refused += 1; continue
            t = pc()
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): T['step'] += pc() - t; N['illegal'] += 1; continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: T['step'] += pc() - t; N['dropped'] += 1; continue
            ck = search._coxeterKeyOrNone(ch); T['step'] += pc() - t; N['step'] += 1
            if mode in 'AG' and ck == base and not tp: FAIL.append((LV[id(alg)], TAINT[id(alg)], str(fingerprint.canonicalKey(alg)), v))
            if mode in 'AG': tab[('J==0' if not J else 'J!=0', 'tp' if tp else 'notp', 'keykept' if ck == base else 'keymoved')] += 1
            if ck == base and not (mode == 'G' and not tp): add(ch, LV[id(alg)] + 1, TAINT[id(alg)] or (mode in 'AGT' and not tp))
    level += 1
tot = pc() - t0
print('n', n, 'cls', cls, 'seeds', len(classes[base]), 'mode', mode, 'expanded', nexp, 'levels', level, 'frontier left', len(frontier), 'total %.1fs' % tot)
print('gate-admitted', N['gate'], 'legal-steps(mutated)', N['step'], 'illegal', N['illegal'], 'dropped', N['dropped'], 'refused', refused)
for k in ('gate', 'step', 'J', 'tp'):
    if N[k]: print('  %-5s total %8.2fs  per call %8.2f ms  (%d calls)' % (k, T[k], 1000 * T[k] / N[k], N[k]))
print('per mutated step (all costs / mutated steps): %.2f ms' % (1000 * tot / max(N['step'], 1)))
for k, c in sorted(tab.items()): print('  ', k, c)
if FAIL:
    print('key-kept tp failures: steps', len(FAIL), 'distinct (parent,v)', len(set((a, b) for _, _, a, b in FAIL)),
          'parent tainted', sum(1 for _, t, _, _ in FAIL if t))
    print('  parent depth histogram', sorted(Counter(l for l, _, _, _ in FAIL).items()))
    print('  parent depth (untainted parents only)', sorted(Counter(l for l, t, _, _ in FAIL if not t).items()))
