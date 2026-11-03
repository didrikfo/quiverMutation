"""Round 111 theorist: attempt the degree-0 loop control: M2-type parent with a degree-0 loop g at i (d = 0), same engine as E-284.
Usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/111/theorist_deg0.py [WMAX]"""
import importlib.util, sys
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)
W = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for wg in (3, 1, 0):
    q = lf.mk({'a': (0, 1, 0), 'b': (1, 0, 0), 'g': (0, 0, 0)}, ['a'], {}, {'a': 1, 'b': 2, 'g': wg})
    print('weight of g =', wg)
    try:
        m, letters, wt = lf.setup(q)
        for L in letters: print('  ', L, 'deg', m.deg(L), 'wt', wt(L), 'd', m.partial({(L,): 1}))
        lf.run('deg0loop', q, W)
    except Exception as e:
        print('  FAILS:', type(e).__name__, str(e)[:200])
