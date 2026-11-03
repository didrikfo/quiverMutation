"""Round 111 theorist, T5/H-023: print all letters of M2 (free degree-1 loop f at i, d f = 0; E-284) with their differentials.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/111/theorist_m2dtheta.py"""
import importlib.util
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)
q = lf.MODELS['M2_loop1_dphi=0'](); m, letters, wt = lf.setup(q)
for L in letters:
    print(L, m.src(L), '->', m.tgt(L), 'deg', m.deg(L), 'wt', wt(L), 'd =', m.partial({(L,): 1}))
