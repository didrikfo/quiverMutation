"""Round 109 theorist: print d(theta) for the theta letters of the loop models of round 106 (M1, M4, M5), to see which Lemma A hypothesis fails.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/109/theorist_dtheta.py"""
import importlib.util, sys
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)
for name in ['M1_loop1_dphi=ba', 'M4_two_alpha', 'M5_loop_two_paths', 'M6_noloop_control']:
    q = lf.MODELS[name](); m, letters, wt = lf.setup(q)
    print('==', name, 'A =', q.A)
    for L in letters:
        if L[0] == 'm' and L[3] is not None:
            print(' theta', L, 'deg', m.deg(L), 'wt', wt(L), ' d =', m.partial({(L,): 1}))
print('== all letters of M1 (kind, a, p, b) src->tgt deg wt, d')
q = lf.MODELS['M1_loop1_dphi=ba'](); m, letters, wt = lf.setup(q)
for L in letters:
    print(' ', L, m.src(L), '->', m.tgt(L), 'deg', m.deg(L), 'wt', wt(L), 'd =', m.partial({(L,): 1}))
