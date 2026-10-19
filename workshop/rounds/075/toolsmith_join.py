"""Round 075 (toolsmith), T10 item 1: J = 0 tilting-only (F/R) join search for the 4 unjoined presentation classes of E-213 (parent and child of the
printed representatives in rounds/074/skeptic_failkids_an_out.txt).  Moves, target ball, meet search = rounds/049/toolsmith_tiltpath.py (exec'd).
  toolsmith_join.py REP DT DC CAP [CACHEDIR]     REP = 0..3 (order of the 'NEW class rep' blocks), DT target-ball depth, DC algebra-ball depth, CAP max keys.
Each side is searched separately (parent, child); prints HIT total (algebra side + target side) or MISS with level sizes and whether CAP stopped it.
Run from the repo root.  A MISS is a bound for these depths/caps, not a proof of inequivalence."""
import sys, re, ast, time, pickle, os
REP, DT, DC, CAP = (int(x) for x in sys.argv[1:5]); CD = sys.argv[5] if len(sys.argv) > 5 else '/tmp/claude-0/ts75'
txt = open('workshop/rounds/074/skeptic_failkids_an_out.txt').read().split('\n')
blocks = []
for i, l in enumerate(txt):
    if l.startswith('NEW class rep'):
        cls, v = re.match(r'NEW class rep: cls (\d) v (\d)', l).groups()
        d = {}
        for j in range(i + 1, i + 5):
            m = re.match(r'\s+(parent|child)\s+arrows (\[.*?\])\s*$', txt[j])
            if m: cur = m.group(1); d[cur] = [ast.literal_eval(m.group(2)), None]; continue
            m = re.match(r'\s+rels (\[.*\]) *$', txt[j])
            if m: d[cur][1] = eval(m.group(1))
        blocks.append((int(cls), int(v), d))
CLS, V0, D = blocks[REP]
sys.argv = ['x', 'none', str(CLS), str(DT), str(DC), str(CAP), 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
import networkx as nx
def build(arrows, rels):
    q = nx.MultiDiGraph()
    for (a, b, k) in arrows: q.add_edge(a, b, key=k)
    return procedure.toPathAlgebra(q, rels)
P, C = build(*D['parent']), build(*D['child'])
print('rep', REP, 'class', CLS, 'v', V0, 'parent key class ok', search._coxeterKeyOrNone(P) == base, 'child', search._coxeterKeyOrNone(C) == base, flush=True)
# sanity: the printed step is gated, tiltingPlus false, and reproduces the child
rels = procedure.relationsFrom(P)
print('step v=%d: gate' % V0, mutation.mutationIsPossibleAtVertex(P, V0), 'perI', perI(P, V0), 'tiltingPlus', tiltingPlus(P.quiver, rels, V0), flush=True)
t0 = time.time()
cache = '%s/tball_c%d_d%d.pkl' % (CD, CLS, DT)
if os.path.exists(cache): tball, tl, tc = pickle.load(open(cache, 'rb'))
else: tball, tl, tc = ball(classes[base], DT, CAP, 'target'); pickle.dump((tball, tl, tc), open(cache, 'wb'))
print('target ball: seeds', len(classes[base]), 'levels', tl, 'keys', len(tball), 'capped', tc, '%.0fs' % (time.time() - t0), flush=True)
for name, alg in (('child', C), ('parent', P)):
    t1 = time.time(); best, lv, capd = reach(alg, tball, DC, CAP)
    print(name, 'HIT total %d (alg side %d + target side %d)' % best if best else 'MISS', 'levels', lv, 'CAPPED' if capd else 'uncapped', '%.0fs' % (time.time() - t1), flush=True)
