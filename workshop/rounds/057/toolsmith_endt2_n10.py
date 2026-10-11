"""Round 057 (toolsmith): End(T) vs next algebra on the 7-step E-169 paths (n = 10, start LNA [0,5,0,4,0,3,3,0]); signed moves s: F|s| on the algebra (s > 0) or on its dual (s < 0), as experimentalist_amerge.replay.
Usage (repo root): .venv/bin/python workshop/rounds/057/toolsmith_endt2_n10.py plan | run LO:HI   (paths are the 5 'full' lists of workshop/rounds/054/experimentalist_amerge_out.txt)"""
import sys, re, time
MODE = sys.argv[1]; SL = sys.argv[2] if len(sys.argv) > 2 else None
exec(open('workshop/rounds/057/toolsmith_endt2.py').read().replace("PK = os.environ.get('PK'", "PK = os.environ.get('PK'").replace("sys.argv = ['x', PK", "_SAVE = sys.argv; sys.argv = ['x', PK"))
sys.argv = _SAVE
START = [0, 5, 0, 4, 0, 3, 3, 0]
paths = [eval(m) for m in re.findall(r'-> full (\[.*?\])', open('workshop/rounds/054/experimentalist_amerge_out.txt').read())]
if MODE == 'plan':
    for i, p in enumerate(paths): print(i, p)
    sys.exit()
lo, hi = (int(x) for x in SL.split(':')); tot = {}
for n_, p in enumerate(paths):
    if not lo <= n_ < hi: continue
    x = nk.LinearNakayamaAlgebra(10, list(START))
    for s in p:
        t0 = time.time(); kind = 'F' if s > 0 else 'R'; v = abs(s)
        (lab, a, v, c), = list(edges_of(x, ['%s%d' % (kind, v)]))
        rr = compare2(a, v, c); vd = rr['verdict'].split(' (')[0]; tot[vd] = tot.get(vd, 0) + 1
        print('n10 path', n_, lab, 'dims', rr['dims'], 'arrows', rr['arrows'], rr['n_arrows'], 'par', rr['par'], 'maxdim', rr['maxdim'], '|', rr['verdict'], '%.1fs' % (time.time() - t0), flush=True)
        x = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
print('TOTAL', tot)
