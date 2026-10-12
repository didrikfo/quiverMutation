"""Round 061 (toolsmith), referee item 3: power of stepTest at n = 10.  For the parent algebra of each of the 35 E-169 edges, apply stepTest at EVERY vertex v' (both a and dual(a)), i.e. the wrong vertex / reversed direction as well as the used step; compare the verdict (Hom(T,T[+-1]) = 0) with tiltingPlus (context only).
Usage (repo root): .venv/bin/python workshop/rounds/061/toolsmith_homn10_power.py"""
import sys, re, time, collections
sys.argv = ['x', 'plan']
src = open('workshop/rounds/061/toolsmith_homn10.py').read().split("if MODE in ('plan', 'n10')")[0]
exec(src)
paths = [eval(m) for m in re.findall(r'-> full (\[.*?\])', open('workshop/rounds/054/experimentalist_amerge_out.txt').read())]
tab = collections.Counter(); rej = []
for i, p in enumerate(paths):
    x = nk.LinearNakayamaAlgebra(10, list(START))
    for j, s in enumerate(p):
        v0 = abs(s); a0 = pathAlgebra.dualPathAlgebra(x) if s < 0 else x
        for dual in (False, True):
            a = pathAlgebra.dualPathAlgebra(a0) if dual else a0
            for v in sorted(a.quiver.nodes):
                used = (v == v0 and not dual)
                try:
                    Vs, res = stepTest(a, v, -1); m1 = int(np.abs(np.array(res[-1])).sum()); p1 = int(np.abs(np.array(res[1])).sum())
                    tp = bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v))
                except Exception as e:
                    tab['error', used] += 1; continue
                ok = (m1 == 0 and p1 == 0); tab['pass' if ok else 'REJECT', 'tp_pass' if tp else 'tp_fail', 'used' if used else 'other'] += 1
                if not ok: rej.append((i, j, s, 'dual' if dual else 'same', v, m1, p1, tp))
        c, r = edge(a0, v0); x = pathAlgebra.dualPathAlgebra(c) if s < 0 else c
for k, v in sorted(tab.items(), key=str): print(k, v)
print('rejections:', len(rej)); [print(r) for r in rej[:40]]
