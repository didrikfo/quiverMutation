"""Round 115 (theorist): HH^0,1,2 of parent and child of E-201's 13 key-kept D = 0 J != 0 steps (n = 7, class 1), from the checkpoint written by
`rounds/069/skeptic_c1_13.py 1 540 23000` (run twice) at /tmp/skeptic_c1_13_cls1.pkl.  Planted corruptions: child with one relation dropped / one term dropped."""
import sys, pickle, copy
import numpy as np
sys.path.insert(0, '.')
exec(open('workshop/rounds/115/theorist_hh.py').read().split("if __name__")[0])
rows = [r for r in pickle.load(open('/tmp/skeptic_c1_13_cls1.pkl', 'rb')) if r['tag'] == 'E147']
def cart(V, E, R):
    A = Alg(V, E, R); ix = {v: i for i, v in enumerate(A.V)}; C = np.zeros((len(V), len(V)), dtype=int)
    for p, s, t in A.basis: C[ix[s], ix[t]] += 1
    return C
print('steps', len(rows))
nyes = 0; ctl = [0, 0, 0]
for i, r in enumerate(rows):
    V = r['V']; pe = [tuple(e) for e in r['pe']]; ce = [tuple(e) for e in r['ce']]
    pr = eval(r['pr']); cr = eval(r['cr'])
    hp = hh(V, pe, pr, top=2); hc = hh(V, ce, cr, top=2)
    same = hp[1] == hc[1]; nyes += same
    print(i, 'v', r['v'], 'par', r['par'], 'dim', hp[0], hc[0], 'HH parent', hp[1], 'HH child', hc[1], 'cochains', hp[2], '->', 'YES' if same else 'NO', flush=True)
    # planted corruptions of the child
    Cc = cart(V, ce, cr)
    for j in range(len(cr)):
        cor = [dict(x) for x in cr]
        if len(cor[j]) > 1:
            k0 = sorted(cor[j])[0]; del cor[j][k0]            # drop one term of a multi-term relation
        else:
            del cor[j]                                         # drop a monomial relation
        try:
            hk = hh(V, ce, cor, top=2); Ck = cart(V, ce, cor)
        except Exception as ex:
            continue
        ctl[0] += 1; eqC = bool((Ck == Cc).all())
        if hk[1] != hc[1]: ctl[1] += 1
        if eqC and hk[1] != hc[1]: ctl[2] += 1
print('YES (HH equal) on', nyes, 'of', len(rows))
print('planted corruptions of children: total', ctl[0], 'HH differs from child', ctl[1], 'of which Cartan equal to child', ctl[2])
