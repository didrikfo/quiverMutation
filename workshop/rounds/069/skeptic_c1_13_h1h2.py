"""Round 069 (skeptic): E-147's 'strict H1/H2' test (shape/strict copied verbatim from rounds/043/skeptic_c2.py) on the first 13 D = 0 steps
of the class-1 walk checkpoint /tmp/skeptic_c1_13_ck1.pkl. Run from repo root: .venv/bin/python workshop/rounds/069/skeptic_c1_13_h1h2.py"""
import sys, pickle
import numpy as np
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]; exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
s = open('workshop/rounds/069/skeptic_c1_13.py').read(); exec(s[s.index("from fractions"):s.index("classes = {}")])
t = open('workshop/rounds/043/skeptic_c2.py').read(); exec(t[t.index("def shape"):t.index("def orbitrel")])
seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
Z = sorted([r for r in byCart.values() if detpoly(r['CA']) == detpoly(r['CB'])], key=lambda r: r['nexp'])
nh = ns = 0
for k, r in enumerate(Z[:13], 1):
    vi = r['V'].index(r['v']); ok, w, o = shape(r['CA'], r['CB'], vi); st = None
    if ok:
        i = [r['V'].index(x) for x in r['J']][0]; st = strict(r['CA'], r['CB'], vi, w, i, o)
    nh += bool(ok); ns += bool(st); print(k, 'H1/H2 shape', ok, 'strict', st, 'J', r['J'])
print('first 13: H1/H2 shape', nh, 'strict', ns)
print('all D=0 Cartan steps in checkpoint:', len(Z), 'shape', sum(1 for r in Z if shape(r['CA'], r['CB'], r['V'].index(r['v']))[0]))
