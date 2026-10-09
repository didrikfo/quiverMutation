"""Round 050 (skeptic): apply the independent tilting test to the recorded E-151/E-154 steps (rounds/049/toolsmith_collect.py pickle):
the 'fail' steps (key-keeping, J != 0 / tiltingPlus false) and the 'ctrl' steps (J = 0).  Usage: skeptic_failsteps.py in.pkl"""
import sys, pickle
import numpy as np
sys.path.insert(0, '.'); PK = sys.argv[1]; sys.argv = ['x']
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import invariants
from collections import Counter
cnt = Counter(); rows = []
for i, r in enumerate(pickle.load(open(PK, 'rb'))):
    a = r['parentObj']; c = r['childObj']; v = r['v']
    V, res = stepTest(a, v, -1)
    m1 = int(np.array(res[-1]).sum()); p1 = int(np.array(res[1]).sum())
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    cart = bool((np.array(res[0]) == C.T).all())
    cnt[(r['kind'], bool(r['J']), bool(r['tp']), m1 == 0 and p1 == 0, cart)] += 1
    rows.append((i, r['kind'], r['depth'], v, bool(r['J']), m1, p1, cart))
for k, c in sorted(cnt.items(), key=str): print('kind %s J!=0 %s tiltingPlus %s indep-tilt %s H==C(child)^T %s : %d' % (k + (c,)))
for x in rows:
    if x[1] == 'fail': print(x)
