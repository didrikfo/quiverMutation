"""T5 round 019: scholar_walk.py (round 014, E-084) with checkpoint/resume. Same walk, same counts, same output
lines; the library is used as it is (E-089 fix is in it), no monkeypatch.
  toolsmith_walk.py n --plan
  toolsmith_walk.py n --class I [--depth D] [--budget-sec S] [--stop-on-reject] [--ckpt FILE]
With --ckpt FILE: if FILE exists the walk resumes from it (n, class, depth must match); when the budget (seconds,
this slice only) is spent the state is written to FILE (atomically) and the script exits 2 printing 'CHECKPOINT';
at the end (closed or depth reached) FILE is kept with done=True and the summary is printed. Exit 0 if closed,
2 if stopped by depth or budget. Checkpoints hold algebras (pickle), keep them under /tmp when big.
Budget is checked before each expansion, so a slice may overrun by one expansion (< 1 s)."""
import argparse, sys, time, os, pickle
from collections import Counter
sys.path.insert(0, '.')
_argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
ap_ = argparse.ArgumentParser()
ap_.add_argument('n', type=int); ap_.add_argument('--plan', action='store_true')
ap_.add_argument('--class', type=int, dest='cls', default=-1); ap_.add_argument('--depth', type=int, default=99)
ap_.add_argument('--stop-on-reject', action='store_true', dest='sor')
ap_.add_argument('--budget-sec', type=float, default=0, dest='budget'); ap_.add_argument('--show', type=int, default=5)
ap_.add_argument('--ckpt', default='')
ap_.add_argument('--max-exp', type=int, default=0, dest='maxexp', help='stop (like a spent budget) when this many expansions are done in total; deterministic')
a = ap_.parse_args()

classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        k = search._coxeterKeyOrNone(alg)
        classes.setdefault(k, []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
if a.plan:
    for i, k in enumerate(order): print(i, len(classes[k]), k)
    sys.exit(0)
starts = classes[order[a.cls]]; base = order[a.cls]
t_start = time.time()
if a.ckpt and os.path.exists(a.ckpt):
    S = pickle.load(open(a.ckpt, 'rb'))
    assert (S['n'], S['cls']) == (a.n, a.cls), 'checkpoint is for another walk'
    print('RESUME from', a.ckpt, 'depth', S['d'] + 1, 'pos', S['pos'], 'expanded', S['expanded'], 'slices', S['slices'], flush=True)
    S['slices'] += 1
else:
    S = dict(n=a.n, cls=a.cls, seen=set(), cur=[], nxt=[], d=0, pos=0, tab=Counter(), rej=[], mism=[], expanded=0,
             nonmono=0, stop=False, done=False, slices=1, elapsed=0.0)
    for i, s in enumerate(starts):
        key = fingerprint.canonicalKey(s)
        if key is not None:
            if key in S['seen']: continue
            S['seen'].add(key)
        S['nxt'].append((s, (i,)))
    S['cur'] = S['nxt']; S['nxt'] = []
seen, tab, rej, mism = S['seen'], S['tab'], S['rej'], S['mism']
def save():
    S['elapsed'] += time.time() - t_start
    tmp = a.ckpt + '.tmp'; pickle.dump(S, open(tmp, 'wb'), protocol=pickle.HIGHEST_PROTOCOL); os.replace(tmp, a.ckpt)
def add(alg, path):
    key = fingerprint.canonicalKey(alg)
    if key is not None:
        if key in seen: return
        seen.add(key)
    S['nxt'].append((alg, path))
t0 = time.time(); budget_hit = False
while not S['done'] and S['d'] < a.depth:
    d = S['d']
    if not S['cur'] and S['pos'] == 0 and d > 0 and not S['nxt']: S['done'] = True; break
    cur = S['cur']
    while S['pos'] < len(cur):
        if (a.budget and time.time() - t0 > a.budget) or (a.maxexp and S['expanded'] >= a.maxexp): budget_hit = True; break
        alg, path = cur[S['pos']]
        S['expanded'] += 1; S['pos'] += 1
        if list(nx.simple_cycles(alg.quiver)): tab['parent-cyclic'] += 1; continue
        rels = procedure.relationsFrom(alg)
        if any(len(r) > 1 for r in rels): S['nonmono'] += 1
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            t = tiltingPlus(alg.quiver, rels, v)
            child = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(child.quiver, r) for r in procedure.relationsFrom(child)):
                tab[('illegal', t)] += 1; continue
            child = reduction.reducePathAlgebra(child)
            ck = search._coxeterKeyOrNone(child)
            guard = ck == base
            tab[('guard' if guard else 'noguard', 'tilt' if t else 'NOTtilt')] += 1
            if t and not guard: mism.append((d, alg.rels, v, path))
            if not t: rej.append((d, alg.rels, v, guard, path))
            if guard: add(child, path + (v,))
    if budget_hit: break
    print('depth', d + 1, 'expanded', S['expanded'], 'next', len(S['nxt']), 'nonmono parents', S['nonmono'], '%.0fs' % (S['elapsed'] + time.time() - t_start), flush=True)
    S['cur'], S['nxt'], S['pos'], S['d'] = S['nxt'], [], 0, d + 1
    if a.sor and rej: S['stop'] = True; break
    if not S['cur']: S['done'] = True
if budget_hit:
    if a.ckpt: save(); print('CHECKPOINT', a.ckpt, 'depth', S['d'] + 1, 'pos', S['pos'], 'of', len(S['cur']), 'expanded', S['expanded'], flush=True)
    else: print('budget spent, no --ckpt: counts partial', flush=True)
closed = S['done'] and not budget_hit and not S['stop']
if a.ckpt and not budget_hit: save()
if True:
    print('n', a.n, 'class', a.cls, 'key', base, 'starts', len(starts), 'distinct algebras', len(seen), 'closed', closed)
    for k, v in sorted(tab.items(), key=str): print(k, v)
    print('REJECTIONS (gate-admitted, tiltingPlus False):', len(rej))
    for r in rej[:a.show]: print(r)
    print('GATE+TILT BUT KEY MOVES:', len(mism))
    for r in mism[:a.show]: print('M', r)
sys.exit(0 if closed else 2)
