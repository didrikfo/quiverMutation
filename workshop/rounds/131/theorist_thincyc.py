"""Round 131, theorist: are the quivers at gate-admitted steps of F-041 balls acyclic?
(The thinness lemma's library reading -- simple-path d, gate, isInIdeal -- is exact on acyclic quivers.)
Reuses round-130 toolsmith_f041calls.py; adds a cyclic counter on isTiltingAtVertex.
usage: theorist_thincyc.py N idx[,idx...]"""
import sys, json, importlib.util, warnings
warnings.simplefilter('ignore')
sys.path.insert(0, '.')
spec = importlib.util.spec_from_file_location('tfc', 'workshop/rounds/130/toolsmith_f041calls.py')
tfc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tfc)
import networkx as nx
from quivermutation import mutation
CYC = [0, 0]
_w = mutation.isTiltingAtVertex
def _c(alg, v):
    CYC[0] += 1; CYC[1] += not nx.is_directed_acyclic_graph(alg.quiver)
    return _w(alg, v)
mutation.isTiltingAtVertex = _c
N = int(sys.argv[1]); _, _, P = tfc.tfr.orbitsAndPairs(N)
for i in map(int, sys.argv[2].split(',')):
    CYC[:] = [0, 0]
    rec = tfc.onePair(N, *P[i], False)
    print(json.dumps(dict(i=i, calls=rec['calls'], d2=rec['d2'], dmax=rec['dmax'], refused=rec['refusedL'] + rec['refusedR'],
                          counted=CYC[0], cyclic=CYC[1])), flush=True)
