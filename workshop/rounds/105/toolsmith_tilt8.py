"""Round 105 (toolsmith), T10 / E-278: independent tilting test (rounds/050 skeptic_tilt.stepTest + det check) on the key-kept J != 0 edges
of E-268 that are reachable by PRINTED paths: A = the E-254 edge (path of rounds/098 toolsmith_j98.py), B = the E-244 edges (6 printed parents).
NOT covered: the E-239 depth-5 B edge (needs the 3304000 depth-5 ball scan, not rebuilt).  n = 9.
For each kept candidate (kind, v) at the parent x: a = x (F) or dual(x) (R); T = T_v(a); verdict hom (Hom(T,T[+-1]) = 0), det (|det H| = |det Cartan a|),
cart (H == Cartan(reduced child before dualising)^T); context: perI (J), gate.  usage: toolsmith_tilt8.py"""
import sys, time
MYA = list(sys.argv); sys.argv = ['x', 'ctl', 'x']
sys.path.insert(0, '.')
src = open('workshop/rounds/087/toolsmith_d6scan.py').read().split("if MODE == 'build'")[0]
exec(compile(src, 's', 'exec'))
_sj = open('workshop/rounds/082/skeptic_join.py').read(); exec(_sj[_sj.index('def child('):_sj.index('def mkroot(')])
import numpy as np
from quivermutation import invariants
B_PATHS = {8063: [('F', 2), ('R', 8), ('R', 7), ('R', 9), ('F', 7), ('R', 9)], 8170: [('F', 2), ('R', 8), ('R', 9), ('F', 2), ('R', 7), ('R', 7)],
           8196: [('F', 2), ('R', 8), ('R', 9), ('F', 4), ('R', 7), ('R', 7)], 8427: [('F', 2), ('R', 9), ('R', 8), ('R', 7), ('R', 9), ('R', 7)],
           11123: [('R', 8), ('R', 9), ('R', 7), ('R', 7), ('F', 2), ('F', 7)], 4606: [('F', 2), ('F', 1), ('R', 8), ('R', 9), ('R', 7), ('R', 7)]}
A_PATH = [('F', 1), ('F', 1), ('F', 2), ('F', 3), ('R', 9), ('R', 7), ('R', 9), ('F', 3)]
def replay(root, mv):
    x = root
    for kind, v in mv: x = child(x, kind, v)
    return x
def run(tag, x):
    for kind, v, st in jcands(x):
        if st != 'kept': continue
        a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        t0 = time.time()
        raw = mutation.quiverMutationAtVertex(a, v); c = reduction.reducePathAlgebra(raw)
        Vs, res = stepTest(a, v, -1)
        H = np.array(res[0], dtype=int)
        hom = not np.array(res[-1]).any() and not np.array(res[1]).any()
        C0 = np.array(invariants.cartanMatrix(a, exact=True).tolist(), dtype=int)
        Cc = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
        print(tag, kind, v, 'J', bool(perI(a, v)), 'gate', bool(mutation.mutationIsPossibleAtVertex(a, v)), '| hom0', hom, 'mass', (int(np.abs(res[-1]).sum()), int(np.abs(res[1]).sum())),
              'det', round(abs(np.linalg.det(H))) == round(abs(np.linalg.det(C0))), 'H==Cartan(child)^T', bool((H == Cc.T).all()), 'H==Cartan(child)', bool((H == Cc).all()),
              '%.0fs' % (time.time() - t0), flush=True)
rootA = nk.LinearNakayamaAlgebra.fromClassName('2400230'); rootB = nk.LinearNakayamaAlgebra.fromClassName('3304000')
run('A_E254', replay(rootA, A_PATH))
for idx, p in B_PATHS.items(): run('B_E244_%d' % idx, replay(rootB, p))
