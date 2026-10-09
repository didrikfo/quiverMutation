"""Round 026 (toolsmith): the class-0 reject walk (scholar_longsquare.py table + experimentalist_rejects.py classifier) with checkpoint/resume.
  toolsmith_rejwalk.py n --plan                          class sizes only (no walk)
  toolsmith_rejwalk.py n --class I [--ckpt FILE] [--budget-hours H | --budget-sec S] [--max-exp N] [--ckpt-every SEC] [--show K]
Same guarded BFS as rounds 023/025 (level by level, canonical-key dedup, children kept iff Coxeter key == base). Per gate-admitted
step it tabulates (J != 0, tiltingPlus, out-degree capped at 2, longSquare NEW (arrow model, handles parallel arrows), longSquareOLD (r023),
mono). Each distinct (parent key, v) with J != 0 is kept once and classified: out-degree >= 2: 'D' / 'D-part k/m' / 'none' (as r025
experimentalist_rejects.py); out-degree 1: 'sq' (new test) / 'sq-parallel-only' (new True, old False = E-108's missed case) / 'nosq'.
Checkpoint: pickle of the whole state, written atomically (tmp + rename) when the budget is spent, on SIGTERM/SIGINT, and every
--ckpt-every seconds (default 900; 0 = only at the end). Resuming reads FILE (n and class must match). Exit codes: 0 walk closed, 2 stopped
(budget, --max-exp, signal; checkpoint written). The budget is per slice and is checked before each expansion (overrun < ~1 s).
Resume == uninterrupted: --max-exp stops at an exact expansion count, so slices are deterministic (test: see round 026 toolsmith.md).
Checkpoints are large at n = 9: keep FILE under /tmp or another untracked dir."""
import argparse, sys, time, os, pickle, signal, importlib.util
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
_s = importlib.util.spec_from_file_location('toolsmith_longsquare', 'workshop/rounds/026/toolsmith_longsquare.py')
LS = importlib.util.module_from_spec(_s); _s.loader.exec_module(LS)
p = argparse.ArgumentParser(); p.add_argument('n', type=int); p.add_argument('--plan', action='store_true')
p.add_argument('--class', type=int, dest='cls', default=0); p.add_argument('--ckpt', default='')
p.add_argument('--budget-hours', type=float, default=0, dest='bh'); p.add_argument('--budget-sec', type=float, default=0, dest='bs')
p.add_argument('--max-exp', type=int, default=0, dest='maxexp'); p.add_argument('--ckpt-every', type=float, default=900, dest='every')
p.add_argument('--show', type=int, default=3); a = p.parse_args()
budget = a.bh * 3600 + a.bs

classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
if a.plan:
    for i, k in enumerate(order): print(i, len(classes[k]), k)
    sys.exit(0)
base = order[a.cls]

def kerdim(alg, v, rels):
    quiver = alg.quiver; outs = ap.arrowsOutOf(quiver, v); tot = 0; mono = False
    for i in quiver.nodes:
        if i == v: continue
        P = ap.allPathsBetween(quiver, i, v)
        if not P: continue
        dimV = len(P) - len(ap.idealBasis(quiver, rels, i, v)); rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): row[(b, kk)] = x
            rows.append(row)
            if not row and not ap.isInIdeal(quiver, rels, ap.combination([q])): mono = True
        tot += dimV - (rank(rows) if rows else 0)
    return tot, mono
def intoArrow(alg, v, b):
    e = b[1]
    return any(len(rel) >= 2 and all(len(q) >= 3 and q[-1] == e and q[-2] == v for q in rel) for rel in alg.rels)
def classify(alg, v, outs, lsNew, lsOld):
    if len(outs) == 1: return ('sq' if lsOld else 'sq-parallel-only') if lsNew else 'nosq'
    k = sum(intoArrow(alg, v, b) for b in outs)
    return 'D' if k == len(outs) else ('D-part%d/%d' % (k, len(outs)) if k else 'none')

if a.ckpt and os.path.exists(a.ckpt):
    S = pickle.load(open(a.ckpt, 'rb')); assert (S['n'], S['cls']) == (a.n, a.cls), 'checkpoint is for another walk'
    S['slices'] += 1
    print('RESUME', a.ckpt, 'level', S['d'], 'pos', S['pos'], 'of', len(S['cur']), 'expanded', S['expanded'], 'seen', len(S['seen']), 'slice', S['slices'], flush=True)
else:
    S = dict(n=a.n, cls=a.cls, seen=set(), cur=[], nxt=[], d=0, pos=0, tab=Counter(), rej={}, expanded=0, done=False, slices=1, elapsed=0.0)
    for s in classes[base]:
        k = fingerprint.canonicalKey(s)
        if k is not None:
            if k in S['seen']: continue
            S['seen'].add(k)
        S['cur'].append(s)
seen, tab, rej = S['seen'], S['tab'], S['rej']
t_start = time.time(); t_last = t_start; stop = []
def save():
    now = time.time(); S['elapsed'] += now - save.t; save.t = now
    tmp = a.ckpt + '.tmp'; pickle.dump(S, open(tmp, 'wb'), protocol=pickle.HIGHEST_PROTOCOL); os.replace(tmp, a.ckpt)
save.t = t_start
for sg in (signal.SIGTERM, signal.SIGINT): signal.signal(sg, lambda *_: stop.append('signal'))
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    S['nxt'].append(alg)
while not S['done']:
    cur = S['cur']
    while S['pos'] < len(cur):
        if stop or (budget and time.time() - t_start > budget) or (a.maxexp and S['expanded'] >= a.maxexp):
            stop.append('budget'); break
        if a.ckpt and a.every and time.time() - t_last > a.every: save(); t_last = time.time()
        alg = cur[S['pos']]; S['expanded'] += 1; S['pos'] += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(alg.quiver, v)
            kd, mono = kerdim(alg, v, rels); t = bool(tiltingPlus(alg.quiver, rels, v))
            lsN, lsO = LS.longSquare(alg, v), LS.longSquareOld(alg, v)
            tab[('J!=0' if kd else 'J=0', 'tiltingPlus', t, 'outdeg', min(len(outs), 2), 'longsq', lsN, 'old', lsO, 'mono', mono)] += 1
            if kd:
                key = (fingerprint.canonicalKey(alg) or str(alg.rels), v)
                if key not in rej: rej[key] = (alg, v, kd, classify(alg, v, outs, lsN, lsO))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
    if stop: break
    print('level', S['d'] + 1, 'expanded', S['expanded'], 'next', len(S['nxt']), 'rejects', len(rej), '%.0fs' % (S['elapsed'] + time.time() - save.t), flush=True)
    S['cur'], S['nxt'], S['pos'], S['d'] = S['nxt'], [], 0, S['d'] + 1
    if not S['cur']: S['done'] = True
if a.ckpt: save()
if stop: print('CHECKPOINT' if a.ckpt else 'STOPPED (no --ckpt: counts partial)', a.ckpt, 'level', S['d'] + 1, 'pos', S['pos'], 'of', len(S['cur']), 'expanded', S['expanded'], flush=True)
print('n', a.n, 'class', a.cls, 'algebras', len(seen), 'expanded', S['expanded'], 'slices', S['slices'], 'closed', S['done'] and not stop, '%.0fs' % S['elapsed'])
for k, x in sorted(tab.items(), key=str): print(k, x)
C = Counter((v[2], v[3]) for v in rej.values())
print('DISTINCT (parent, v) with J != 0:', len(rej))
for k, x in sorted(C.items(), key=str): print('  dim J', k[0], 'class', k[1], x)
for key, (alg, v, kd, cl) in list(rej.items())[:a.show]: print('  e.g. v', v, 'dimJ', kd, cl, 'arrows', sorted(alg.quiver.edges(keys=True)), 'rels', alg.rels)
sys.exit(0 if S['done'] and not stop else 2)
