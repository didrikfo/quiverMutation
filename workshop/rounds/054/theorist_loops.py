"""Round 054 (theorist), response to referee item 2: along every printed path edge (E-157 child + parent paths, E-160 three paths, E-163 path),
check that the algebra a (x for an F move, dual(x) for an R move) has NO LOOP at the mutated vertex v (= all targets of arrowsOutOf(a,v) differ from v),
and also that its quiver has no loop and no oriented cycle at all (H1).  Also checks the final (meeting) algebra of each side is loopless/acyclic.
Usage (repo root): theorist_loops.py in.pkl cls     (in.pkl from rounds/050/toolsmith_collect.py 7 cls 20000 in.pkl 100)
Edges: E-157 from rounds/049/toolsmith_paths_logs.txt (sections paths_c<cls>, parents_c<cls>); E-160 and E-163 hard-coded from the entries."""
import sys, re, pickle
import networkx as nx
A_ = sys.argv; PK, CLS = A_[1], int(A_[2])
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
from quivermutation import arrowPaths as ap
def step(x, kind, v):
    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    V = sorted(x.quiver.nodes)
    c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
    if sorted(c.vertices()) != V: return None
    return (c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)), a
def loopinfo(a, v):
    q = a.quiver
    outs = ap.arrowsOutOf(q, v)
    loop_v = sum(1 for b in outs if b[1] == v)
    loops_any = sum(1 for _ in nx.selfloop_edges(q))
    cyc = any(True for _ in nx.simple_cycles(q))
    return loop_v, loops_any, cyc, len(outs)
tot = dict(paths=0, edges=0, loop_at_v=0, loops_any=0, cyclic=0, noout=0, nokind=0, end_bad=0, nostep=0)
bad = []
def run(label, x0, mv):
    x = x0
    for s_ in mv:
        kind, v = s_[0], int(s_[1:])
        lv, la, cy, no = loopinfo(x if kind == 'F' else pathAlgebra.dualPathAlgebra(x), v)
        tot['edges'] += 1; tot['loop_at_v'] += (lv > 0); tot['loops_any'] += (la > 0); tot['cyclic'] += cy; tot['noout'] += (no == 0)
        if lv or la or cy: bad.append((label, s_, lv, la, cy))
        r = step(x, kind, v)
        if r is None: tot['nostep'] += 1; bad.append((label, s_, 'NOSTEP')); return None
        x = r[0]
    la = sum(1 for _ in nx.selfloop_edges(x.quiver)); cy = any(True for _ in nx.simple_cycles(x.quiver))
    if la or cy: tot['end_bad'] += 1; bad.append((label, 'END', la, cy))
    return x
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
seedkeys = classes[base]
SEC = None
for line in open('workshop/rounds/049/toolsmith_paths_logs.txt'):
    if line.startswith('=='): SEC = line.strip('= \n'); continue
    if SEC is None or not SEC.endswith('_c%d' % CLS): continue
    m = re.match(r'(paths|parents) (\d+) key (\w+) parentdepth (\d+) replay ok total (\d+) \| child->meet: (.*?) \| LNA/dual #(\d+) ->meet: (.*?) \d+s', line)
    if not m: continue
    which, i, key, pd, total, cs, lna, ts = m.groups(); i = int(i); lna = int(lna)
    r = recs[i]; x0 = r['childObj'] if which == 'paths' else r['parentObj']
    run('E155 %s %d child side' % (which, i), x0, cs.split()); run('E155 %s %d LNA side' % (which, i), seedkeys[lna], ts.split()); tot['paths'] += 1
extra = {1: [('E158 c1 child 5', 5, 'F1 F3 F5 F4 R7 F3 R1', 9, 'F7 R2 R1 F2 F7'), ('E158 c1 child 12', 12, 'R7 R4 R7 R1', 6, 'F1 F5 F4 R7 F3'),
             ('E161 c1 child 14', 14, 'F1 F3 F1 F5 F4 R7 R1', 9, 'F7 R2 R1 F2 F7 R3')],
         2: [('E158 c2 child 6', 6, 'F1 F7 R3 R5 R5 F4 R6', 0, 'F1 F1 F2 R6 R5')]}[CLS]
for lab, ci, cs, lna, ts in extra:
    run(lab + ' child side', recs[ci]['childObj'], cs.split()); run(lab + ' LNA side', seedkeys[lna], ts.split()); tot['paths'] += 1
print('SUMMARY cls', CLS, tot); print('BAD', bad)
