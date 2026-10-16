"""Round 069 (skeptic): prefix counts of the class-1 walk from the checkpoint written by skeptic_c1_13.py (/tmp/skeptic_c1_13_ck1.pkl).
Run from repo root: .venv/bin/python workshop/rounds/069/skeptic_c1_13_prefix.py"""
import sys, pickle
sys.path.insert(0, '.'); sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]; exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
s = open('workshop/rounds/069/skeptic_c1_13.py').read(); exec(s[s.index("from fractions"):s.index("classes = {}")])
seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
Z = sorted([r for r in byCart.values() if detpoly(r['CA']) == detpoly(r['CB'])], key=lambda r: r['nexp'])
for k, r in enumerate(Z[:13], 1): print(k, 'first seen at expansion', r['nexp'], 'depth', r['depth'], 'v', r['v'], 'parallel', r['par'])
for N in (20000, 22449, 22450, 24032, 24033):
    print(N, 'distinct J!=0 Cartan steps', sum(r['nexp'] <= N for r in byCart.values()), 'D0', sum(r['nexp'] <= N for r in Z),
          'key-kept (parent,v) pairs', sum(1 for r in byKey.values() if r['nexp'] <= N and r['kept']))
