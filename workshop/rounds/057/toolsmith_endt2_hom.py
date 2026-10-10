"""Round 057 (toolsmith), response to referee item 1: Hom(T,T[-1]), Hom(T,T[1]) at every failing step of class CLS with the skeptic's stepTest (rounds/050/skeptic_tilt.py).  CLS=2 .venv/bin/python workshop/rounds/057/toolsmith_endt2_hom.py"""
import pickle
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
_st = open('workshop/rounds/050/skeptic_tilt.py').read(); P = 2147483647; exec(_st[_st.index('def rank'):])
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']; nz = 0
for n_, r in enumerate(recs):
    Vs, res = stepTest(r['parentObj'], r['v'], -1); m1 = sum(map(sum, res[-1])); p1 = sum(map(sum, res[1])); nz += m1 != 0
    print('c%d' % CLS, n_, 'v', r['v'], 'Hom(T,T[-1])', m1, 'Hom(T,T[1])', p1)
print('TOTAL c%d: %d of %d have Hom(T,T[-1]) != 0' % (CLS, nz, len(recs)))
