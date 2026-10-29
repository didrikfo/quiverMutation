"""Round 098 (toolsmith): T10, 2400230 (F/R ball, depth 8) vs 3304000 (depth 6) at n = 9 with the key-kept J != 0 (perI) edges added as moves on top of F/R and lnaMoves.
Balls: state pickles of workshop/rounds/087/toolsmith_d6scan.py build (DEPTH=8 ... build 2400230 A8.pkl ; DEPTH=6 ... build 3304000 B6.pkl; exit 2 = rerun).
usage (repository root):
  toolsmith_j98.py --plan
  toolsmith_j98.py jedges A8.pkl B6.pkl OUT.pkl   find the key-kept J != 0 steps: A = replay E-254's parent path (F1 F1 F2 F3 R9 R7 R9 F3) and read the kept candidate;
        B = the six E-244 parents (paths from the 087 log) + a scan of all depth-5 B parents (E-239); each child + its radius-1 J = 0 ball is saved
  toolsmith_j98.py meet A8.pkl B6.pkl OUT.pkl     iso meet with the J != 0 children and their radius-1 balls added to their own sides, plus lnaMoves component check of every extra LNA node
J != 0 edges are used as UNDIRECTED edges (the step parent -> child is a derived equivalence only under the J = 0 premise's complement; this is the calibration's assumption, not a theorem)."""
import sys, time, pickle, collections
MYA = list(sys.argv); MODE = MYA[1] if len(MYA) > 1 else '--plan'
if MODE == '--plan':
    print("""sizing (round 098): this IS the round-093 job plus the J != 0 children. Builds A8 (20 210 nodes, ~5 min, 2 commands) and B6 (15 401 nodes, ~11 min, 2 commands) are not committed, so they are rebuilt (parallel, 2 of 4 cores).
 jedges: replay 8 moves + jcands (seconds); 6 B parents replayed (seconds); depth-5 scan of 2 958 B parents ~ 2 min; children + radius-1 balls ~ 1 min.   meet: WL signature of ~35 600 graphs + iso on buckets ~ 6 min (E-254 form).
 every command < 10 min, total about 25 min wall, nothing overnight.""")
    sys.exit(0)
sys.argv = ['x', 'ctl', 'x']
src = open('workshop/rounds/087/toolsmith_d6scan.py').read().split("if MODE == 'build'")[0]
sys.path.insert(0, '.'); exec(compile(src, 's', 'exec'))
MODE = MYA[1]
from quivermutation import lnaMoves, lines
_sj = open('workshop/rounds/082/skeptic_join.py').read(); exec(_sj[_sj.index('def child('):_sj.index('def mkroot(')])  # child, struct, iso, sig
T0 = time.time()
def replay(root, mv):
    x = root
    for kind, v in mv: x = child(x, kind, v)
    return x
def lnaName(x):
    q = x.quiver
    if q.number_of_edges() != len(q.nodes) - 1 or any(q.out_degree(v) > 1 or q.in_degree(v) > 1 for v in q.nodes): return None
    try: _, rl, _ = lnaMoves.standardise(x)
    except Exception: return None
    return None if rl is None else lines.className(rl)
B_PATHS = {8063: [('F', 2), ('R', 8), ('R', 7), ('R', 9), ('F', 7), ('R', 9)], 8170: [('F', 2), ('R', 8), ('R', 9), ('F', 2), ('R', 7), ('R', 7)],
           8196: [('F', 2), ('R', 8), ('R', 9), ('F', 4), ('R', 7), ('R', 7)], 8427: [('F', 2), ('R', 9), ('R', 8), ('R', 7), ('R', 9), ('R', 7)],
           11123: [('R', 8), ('R', 9), ('R', 7), ('R', 7), ('F', 2), ('F', 7)], 4606: [('F', 2), ('F', 1), ('R', 8), ('R', 9), ('R', 7), ('R', 7)]}
A_PATH = [('F', 1), ('F', 1), ('F', 2), ('F', 3), ('R', 9), ('R', 7), ('R', 9), ('F', 3)]
def keptChildren(par):
    return [(k, v, child(par, k, v)) for k, v, st in jcands(par) if st == 'kept']
if MODE == 'jedges':
    SA, SB = pickle.load(open(MYA[2], 'rb')), pickle.load(open(MYA[3], 'rb'))
    rootA = SA['seen'][next(iter(SA['seen']))][0]; rootB = SB['seen'][next(iter(SB['seen']))][0]
    out = {'A': [], 'B': []}
    pa = replay(rootA, A_PATH); print('A parent key in ball:', fingerprint.canonicalKey(pa) in SA['seen'], 'depth', SA['seen'].get(fingerprint.canonicalKey(pa), (0, 0, 0, None))[3], flush=True)
    for k, v, c in keptChildren(pa): out['A'].append(('E254', k, v, c)); print('A kept J != 0 step', k, v, 'LNA', lnaName(c), flush=True)
    for idx, p in B_PATHS.items():
        pb = replay(rootB, p)
        for k, v, c in keptChildren(pb): out['B'].append(('E244_%d' % idx, k, v, c)); print('B', idx, 'kept', k, v, 'LNA', lnaName(c), flush=True)
    seen6 = {fingerprint.canonicalKey(c) for _, _, _, c in out['B']}
    n5 = 0
    for key, (x, pk, mv, d) in list(SB['seen'].items()):
        if d != 5: continue
        n5 += 1
        for k, v, c in keptChildren(x):
            out['B'].append(('d5', k, v, c)); print('B depth-5 kept', k, v, 'LNA', lnaName(c), flush=True)
    print('depth-5 B parents scanned', n5, 'time %.0fs' % (time.time() - T0), flush=True)
    pickle.dump(out, open(MYA[4], 'wb'))
elif MODE == 'meet':
    SA, SB = pickle.load(open(MYA[2], 'rb')), pickle.load(open(MYA[3], 'rb')); J = pickle.load(open(MYA[4], 'rb'))
    La, Lb = dict(SA['seen']), dict(SB['seen']); extraNodes = []
    for side, L, dd in (('A', La, 9), ('B', Lb, 7)):
        for i, (tag, k, v, c) in enumerate(J[side]):
            L[('x', tag, i)] = (c, None, None, dd); extraNodes.append((side, tag, c))
            for k2, v2, c2 in moves(c): L[('x', tag, i, k2, v2)] = (c2, None, None, dd + 1); extraNodes.append((side, tag + '+1', c2))
    print('extras: A', sum(1 for e in extraNodes if e[0] == 'A'), 'B', sum(1 for e in extraNodes if e[0] == 'B'), flush=True)
    # lnaMoves component check of the extra LNA nodes against the two roots' components (undirected, all rows of length 9)
    par = {}
    def f(a):
        while par.setdefault(a, a) != a: par[a] = par[par[a]]; a = par[a]
        return a
    for r in nk.allRelationLengths(9):
        u = ''.join(map(str, r))
        for w in lnaMoves.movesByRule(9, [int(c) for c in u]): par[f(u)] = f(w)
    fa, fb = f('2400230'), f('3304000'); print('lnaMoves components of the roots equal?', fa == fb, flush=True)
    cnt = collections.Counter()
    for side, tag, c in extraNodes:
        nm = lnaName(c)
        if nm is not None: cnt[(side, 'LNA', 'compA' if f(nm) == fa else 'compB' if f(nm) == fb else 'other')] += 1
    print('extra LNA nodes by side and lnaMoves component:', dict(cnt), flush=True)
    buck = {}
    for k, v in Lb.items(): buck.setdefault(sig(v[0]), []).append(k)
    best = None; bestx = None; ncmp = nb = 0
    for ka, va in La.items():
        for kb in buck.get(sig(va[0]), []):
            nb += 1
            if best is not None and va[3] + Lb[kb][3] >= best[0]: continue
            ncmp += 1
            if iso(va[0], Lb[kb][0]): best = (va[3] + Lb[kb][3], ka, kb)
    print('ball sizes', len(La), len(Lb), 'WL-bucket pairs', nb, 'iso tests', ncmp, 'meet', best if best is None else (best[0], str(best[1])[:60], str(best[2])[:60]), 'time %.0fs' % (time.time() - T0))
elif MODE == 'control':   # positive control: the J != 0 children on the B side meet the plain B ball (E-244 joins them at total 6-8), and the A-side child meets the plain A ball (E-254: total 9)
    SA, SB = pickle.load(open(MYA[2], 'rb')), pickle.load(open(MYA[3], 'rb')); J = pickle.load(open(MYA[4], 'rb'))
    for side, S in (('A', SA), ('B', SB)):
        buck = {}
        for k, v in S['seen'].items(): buck.setdefault(sig(v[0]), []).append(k)
        for tag, k, v, c in J[side]:
            hit = [(S['seen'][kb][3], kb) for kb in buck.get(sig(c), []) if iso(c, S['seen'][kb][0])]
            print(side, tag, 'child iso to a node of its own plain ball:', min(hit)[0] if hit else None, flush=True)
