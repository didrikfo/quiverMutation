"""Round 063 (toolsmith), T10 item 1 (E-186 open item): the E-180 dimension assertion (compare2 -> dimcheck: dim K Q/I' from crels == dim End(T) pair by pair)
on the printed n = 7 path edges (E-157/E-160/E-163 plans of classes 1 and 2) and the 35 n = 10 E-169 edges.
Usage (repo root; n = 7 needs /tmp/tsm/c<CLS>.pkl, see rounds/057/toolsmith_endt2_run.py header):
  CLS=1 timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges.py n7 LO:HI     (path indices of the plan list; 'plan' lists them)
  timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges.py n10
Per edge prints dimcheck (None = assertion holds), the dimension totals, and compare2's verdict; an exception is reported as UNEVALUABLE with its text."""
import sys, re, time, collections
MODE0 = sys.argv[1]; SL0 = sys.argv[2] if len(sys.argv) > 2 else None
if MODE0 == 'n10':
    _s = sys.argv; src = open('workshop/rounds/057/toolsmith_endt2_n10.py').read().split("START = ")[0]
    sys.argv = [_s[0], 'plan']; exec(src.replace("MODE = sys.argv[1]; SL = sys.argv[2] if len(sys.argv) > 2 else None", "MODE = 'x'; SL = None"))
else:
    _s = sys.argv; sys.argv = [_s[0], 'plan']
    src = open('workshop/rounds/057/toolsmith_endt2_run.py').read().split("P = plan()")[0]; exec(src)
    sys.argv = _s
MODE, SL = MODE0, SL0
def check(a, v, c):
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); A, d, carrows, crels = child_info(c)
    msg = dimcheck(carrows, crels, d); qd = quotient_dims(carrows, crels)
    sameD = all(dims[(i, j)] == d[(i, j)] for i in TE.V for j in TE.V)
    return msg, sum(qd.values()), sum(dims.values()), sum(d.values()), sameD, len(crels)
tot = collections.Counter(); bad = []
def one(label, a, v, c):
    t0 = time.time()
    try:
        msg, q, e, dc, sd, nr = check(a, v, c)
        k = 'holds' if msg is None else 'FAILS'
        if msg is not None or not sd: bad.append((label, msg, sd))
        tot[k + ('' if sd else ' (End dims != child dims)')] += 1
        print(label, 'dimcheck', msg, '| dim Q/I', q, 'dim End(T)', e, 'dim child', dc, 'nrels', nr, '%.1fs' % (time.time() - t0), flush=True)
    except Exception as ex:
        tot['UNEVALUABLE'] += 1; bad.append((label, repr(ex)[:200], None)); print(label, 'UNEVALUABLE', repr(ex)[:200], flush=True)
if MODE0 == 'n10':
    START = [0, 5, 0, 4, 0, 3, 3, 0]
    paths = [eval(m) for m in re.findall(r'-> full (\[.*?\])', open('workshop/rounds/054/experimentalist_amerge_out.txt').read())]
    for n_, p in enumerate(paths):
        x = nk.LinearNakayamaAlgebra(10, list(START))
        for s in p:
            kind = 'F' if s > 0 else 'R'; (lab, a, v, c), = list(edges_of(x, ['%s%d' % (kind, abs(s))]))
            one('n10 path %d %s' % (n_, lab), a, v, c); x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
else:
    P = plan()
    if MODE == 'plan':
        for n_, (kind, i, cm, k, lm) in enumerate(P): print(n_, kind, 'idx', i, 'edges', len(mv(cm)) + len(mv(lm)))
        print('paths', len(P), 'edges', sum(len(mv(p[2])) + len(mv(p[4])) for p in P)); sys.exit()
    lo, hi = (int(x) for x in SL.split(':'))
    recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
    for n_, (kind, i, cm, k, lm) in enumerate(P):
        if not lo <= n_ < hi: continue
        r = recs[i]; start = r['parentObj'] if kind == 'E155 parent' else r['childObj']
        for side, st, ms in (('child', start, mv(cm)), ('lna%d' % k, classes[base][k], mv(lm))):
            for lab, a, v, c in edges_of(st, ms): one('c%d path %d %s %s %s' % (CLS, n_, kind, side, lab), a, v, c)
print('TOTAL', dict(tot)); print('BAD', bad)
