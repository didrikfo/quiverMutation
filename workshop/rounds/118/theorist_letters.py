"""Round 118 (theorist): print all letters of kM (kind,a,p,b), src->tgt, degree, weight, d for the loop models M4, M4' (basis change a = a1 - a2), M5 and
M5' (automorphism a -> a + m e).  Usage: .venv/bin/python -W ignore workshop/rounds/118/theorist_letters.py"""
import importlib.util
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)
lf.MODELS["M4'"] = lambda: lf.mk({'a1': (0, 1, 0), 'a': (0, 1, 0), 'b': (1, 0, 0), 'f': (0, 0, 1)}, ['a1', 'a'], {'f': {('a', 'b'): 1}}, {'a1': 1, 'a': 1, 'b': 2, 'f': 3})
lf.MODELS["M5'"] = lambda: lf.mk({'a': (0, 1, 0), 'b': (1, 0, 0), 'f': (0, 0, 1), 'm': (0, 2, 0), 'e': (2, 1, 0)}, ['a', 'm'], {'f': {('a', 'b'): 1}}, {'a': 2, 'b': 2, 'f': 4, 'm': 1, 'e': 1})
import sys
for name in sys.argv[1:] or ['M4_two_alpha', "M4'", 'M5_loop_two_paths', "M5'"]:
    q = lf.MODELS[name](); m, letters, wt = lf.setup(q)
    print('==', name)
    for L in letters:
        print(' ', L, m.src(L), '->', m.tgt(L), 'deg', m.deg(L), 'wt', wt(L), ' d =', m.partial({(L,): 1}))
