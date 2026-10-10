"""Round 057 (toolsmith): power control of symcheck2 on the 8 parallel-arrow failing c1 steps.  Perturb each relation of c: (B) a multi-term relation replaced by its first term; (A) first coefficient doubled; (C) relation dropped
(power against too many solutions: dropping gives iso trivially, shown as a sanity check only).  Expected: B rejected always; A rejected when the scalars cannot absorb it.  CLS=1 .venv/bin/python workshop/rounds/057/toolsmith_endt2_ctrl.py"""
import pickle
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
tot = dict(A_no=0, A_iso=0, B_no=0, B_iso=0, n=0)
for n_, r in enumerate(recs):
    a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
    carrows = ap.arrowsOf(c.quiver)
    cc = {}
    for (s, h, k) in carrows: cc[(s, h)] = cc.get((s, h), 0) + 1
    if max(cc.values()) < 2: continue
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); crels = procedure.relationsFrom(c); out = []
    for k_, rl in enumerate(crels):
        if len(rl) < 2: continue
        k0 = next(iter(rl)); ra = dict(rl); ra[k0] = ra[k0] * 2; pa = list(crels); pa[k_] = ra; pb = list(crels); pb[k_] = {k0: rl[k0]}
        va = symcheck2(TE, arrows, 'same', carrows, pa); vb = symcheck2(TE, arrows, 'same', carrows, pb)
        tot['A_no' if va.startswith('NO') else 'A_iso'] += 1; tot['B_no' if vb.startswith('NO') else 'B_iso'] += 1; out.append((k_, va[:2], vb[:2]))
    tot['n'] += 1; print('step', n_, 'v', v, 'multi-term relations (idx, A, B):', out, flush=True)
print('CONTROL', tot)
