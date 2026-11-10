import sys
from quivermutation import lnaMoves, edgeMoves, doubleMutation, piecewiseHereditary as pwh, coxeterTables as ct
P15 = "0003000000000 2003000000000 0203000000000 0203000000002 0334000000002 0334000020000 3350000020000 3350000020002 3600000020002 3600000202000 4660000202000 3560003002000 4555003002000 3450040002000 4444040002000 3340500002000 0333500002000 0236000002000 0246600002000 0235600030000 0245550030000 0234500400000 0244440400000 0204033400000 0203335000000 2003335000000 0003335000000".split()
P13 = "00300000000 00300000002 00300002000 20300002000 33400002000 44440002000 50403002000 55040002000 50500020000 45500300000 35004000000 33404000000 20350000000 00350000000".split()
def fam(n, a, b):
    out = []
    for nm, f in (('rule', lambda r: [tuple(x) for x in lnaMoves.rewritesOf(n, list(r))]),
                  ('edge', lambda r: [tuple(x) for x, _ in edgeMoves.rewritesOf(n, r)]),
                  ('dbl', lambda r: [tuple(x) for x, _ in doubleMutation.rewritesOf(n, r)])):
        if b in f(a): out.append(nm+'>')
        if a in f(b): out.append(nm+'<')
    return out
def img(n, r):
    L, R = pwh.removeVertex(n, list(r), 1)
    return ct.lnaCoxeterKey(L, R)
for n, P in ((15, P15), (13, P13)):
    rows = [tuple(map(int, s)) for s in P]
    keys = {ct.lnaCoxeterKey(n, list(r)) for r in rows}
    fams = [fam(n, a, b) for a, b in zip(rows, rows[1:])]
    print(n, 'source keys distinct:', len(keys), 'steps:', len(fams), 'unjustified:', sum(1 for f in fams if not f))
    print('  families:', [f[0] for f in fams])
    print('  ends image keys equal?', img(n, rows[0]) == img(n, rows[-1]))
    print('  ', img(n, rows[0])); print('  ', img(n, rows[-1]))
