"""Round 061 (toolsmith), T10 item 1: Hom(T,T[m]), m = -1,0,1, and Cartan(End T) vs Cartan(child), edge by edge, with the independent test of
workshop/rounds/050/skeptic_tilt.py (stepTest; no tiltingPlus/perI/gate for the verdict) on
  (a) the 5 seven-step n = 10 E-169 paths 05040330 -> 33460000 (paths of workshop/rounds/054/experimentalist_amerge_out.txt);
  (b) positive control: J != 0 failing steps of n = 7 (pickle from workshop/rounds/050/toolsmith_collect.py 7 <cls> 20000 PK 20).
Usage (repo root):  toolsmith_homn10.py plan | n10 [LO:HI] | ctrl PK
tiltingPlus (rounds/001/scholar_h015.py) is printed as context only."""
import sys, re, time, copy, pickle
import numpy as np
sys.path.insert(0, '.')
from quivermutation import nakayama as nk, mutation, procedure, reduction, pathAlgebra, invariants
_st = open('workshop/rounds/050/skeptic_tilt.py').read(); exec(_st[_st.index('from fractions'):])
ns = {}; exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'), ns); tiltingPlus = ns['tiltingPlus']
MODE = sys.argv[1]
START = [0, 5, 0, 4, 0, 3, 3, 0]
def edge(a, v):
    """a = algebra whose vertex v is mutated; returns child c, and the dict of tests."""
    c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(a), v))
    Vs, res = stepTest(a, v, -1)
    m1 = int(np.abs(np.array(res[-1])).sum()); p1 = int(np.abs(np.array(res[1])).sum())
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int); H = np.array(res[0], dtype=int)
    J = bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v))
    return c, dict(m1=m1, p1=p1, cart=bool((H == C.T).all() and sorted(c.vertices()) == Vs), tp=J, ok=(m1 == 0 and p1 == 0))
if MODE in ('plan', 'n10'):
    paths = [eval(m) for m in re.findall(r'-> full (\[.*?\])', open('workshop/rounds/054/experimentalist_amerge_out.txt').read())]
    if MODE == 'plan':
        print(len(paths), 'paths,', sum(map(len, paths)), 'edges'); [print(i, p) for i, p in enumerate(paths)]; sys.exit()
    lo, hi = (int(x) for x in sys.argv[2].split(':')) if len(sys.argv) > 2 else (0, len(paths))
    tot = dict(edges=0, hom0=0, cart=0, both=0, tiltingPlusPass=0)
    for i, p in enumerate(paths[lo:hi], lo):
        x = nk.LinearNakayamaAlgebra(10, list(START))
        for s in p:
            v = abs(s); t0 = time.time(); a = pathAlgebra.dualPathAlgebra(x) if s < 0 else x
            c, r = edge(a, v); x = pathAlgebra.dualPathAlgebra(c) if s < 0 else c
            tot['edges'] += 1; tot['hom0'] += r['m1'] == 0 and r['p1'] == 0; tot['cart'] += r['cart']; tot['both'] += r['ok'] and r['cart']; tot['tiltingPlusPass'] += r['tp']
            print('path', i, 'step', s, 'Hom(T,T[-1])', r['m1'], 'Hom(T,T[1])', r['p1'], 'cartan', r['cart'], 'tiltingPlus', r['tp'], '%.1fs' % (time.time() - t0), flush=True)
    print('TOTAL', tot)
else:
    recs = [r for r in pickle.load(open(sys.argv[2], 'rb')) if r['kind'] == 'fail']; tot = dict(n=0, rejected=0, cart=0, tp=0)
    for k, r in enumerate(recs):
        t0 = time.time(); c, e = edge(r['parentObj'], r['v']); tot['n'] += 1; tot['rejected'] += not e['ok']; tot['cart'] += e['cart']; tot['tp'] += e['tp']
        print('fail', k, 'v', r['v'], 'Hom(T,T[-1])', e['m1'], 'Hom(T,T[1])', e['p1'], 'cartan', e['cart'], 'tiltingPlus', e['tp'], '%.1fs' % (time.time() - t0), flush=True)
    print('TOTAL', tot)
