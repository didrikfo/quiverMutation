"""Round 054: join stage of experimentalist_amerge.py.  Usage (repo root):
.venv/bin/python workshop/rounds/054/experimentalist_amerge_join.py F.pkl B.pkl [maxHits]"""
import sys, pickle, importlib.util
sys.path.insert(0, '.')
spec = importlib.util.spec_from_file_location('am', 'workshop/rounds/054/experimentalist_amerge.py'); am = importlib.util.module_from_spec(spec); spec.loader.exec_module(am)
from quivermutation import search, lnaMoves
F = pickle.load(open(sys.argv[1], 'rb')); B = pickle.load(open(sys.argv[2], 'rb'))
hits = am.join(F, B)
print('fwd', len(F), 'back targets', len(B), 'iso meetings', len(hits), 'shortest total', hits[0][0] if hits else None)
if hits:
    tp = am.loadTP(); seen = 0
    for tot, m, pf, pb, ib, kx, ky, phi in hits[:int(sys.argv[3]) if len(sys.argv) > 3 else 3]:
        path = pf + ib
        alg, out = am.replay(path, tp)
        end = lnaMoves.asRelLengths(alg, 10)
        print('target', ''.join(map(str, m)), 'fwd', pf, 'back', pb, '-> full', path)
        print('  edges', len(out), 'gate', sum(o[1] for o in out), 'J0', sum(o[2] for o in out), 'key', sum(o[3] for o in out), 'end relLengths', end)
