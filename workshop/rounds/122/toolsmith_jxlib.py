"""Round 122 (toolsmith): does the LIBRARY's own guarded walk (search.linesReachedFrom: right mutations from the LNA and from its relation dual,
gate + Coxeter-key guard) from 003033030 (n = 11, 13-LNA SNF(S)-class (..,2,4)) report an LNA of the 59-LNA class (..,2,2)?  Prints every line found
with its SNF(C + C^T) tail.  usage (repo root): DEPTH=5 timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxlib.py [ROOT]"""
import sys, os, time
sys.path.insert(0, '.')
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import search, nakayama as nk, coxeterTables as ct
ROOT = sys.argv[1] if len(sys.argv) > 1 else '003033030'; D = int(os.environ.get('DEPTH', '5')); N = len(ROOT) + 2
def snf(rl):
    C = sp.Matrix(ct.lnaCartanMatrix(N, rl)); S = smith_normal_form(C + C.T, domain=sp.ZZ); return tuple(sorted(abs(int(S[i, i])) for i in range(N)))[-3:]
t0 = time.time(); found = search.linesReachedFrom(nk.LinearNakayamaAlgebra.fromClassName(ROOT), D)
print('root', ROOT, 'SNF tail', snf(tuple(int(c) for c in ROOT)), 'depth', D, '%.0fs' % (time.time() - t0))
for item in (found.items() if isinstance(found, dict) else found):
    print('  ', item if not isinstance(item, tuple) else item[0], flush=True)
