"""Round 137 (toolsmith), review item 2 pre-check: rebuild the starts of E-316's 19 key-kept J != 0 steps (toolsmith_t10walk.e316Start),
and check: start key == root key, gate admits the step's vertex, `isTiltingAtVertex` False (so the option refuses it at step 1), child
key kept, child class (SNF(C + C^T) and signature) equal to the start's.  usage: .venv/bin/python -W ignore workshop/rounds/137/toolsmith_e316check.py"""
import sys, time
sys.path.insert(0, 'workshop/rounds/137'); sys.argv = ['x']
import toolsmith_t10walk as T
from quivermutation import search, nakayama as nk
t0 = time.time()
for i, (rec, root, path, step) in enumerate(T.E316):
    S = T.e316Start(i); v = int(step[1:]); base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName(root))
    legal, c = T.childOf(S, v)
    print('%2d %-12s start key = root key %s  gate %s  tilting %s  legal %s  child key kept %s  child in-class %s'
          % (i, rec, search._coxeterKeyOrNone(S) == base, T._gate(S, v), T._tilt(S, v), legal,
             legal and search._coxeterKeyOrNone(c) == base, legal and T.inv(c) == T.inv(S)), flush=True)
print('time %.0fs' % (time.time() - t0))
