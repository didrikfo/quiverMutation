"""Round 063 (toolsmith): power control for the dimension assertion on path edges: drop each single relation of crels on every edge of class CLS paths
and count dimcheck rejections (a drop is 'redundant' iff dim K Q/I'' equals dim End(T) -> accept is then correct).  Usage: CLS=2 timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges_ctrl.py"""
import sys, collections
sys.argv = ['x', 'n7', '0:0']
exec(open('workshop/rounds/063/toolsmith_dimedges.py').read().split("tot = collections.Counter()")[0])
P = plan(); recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']; T = collections.Counter()
for n_, (kind, i, cm, k, lm) in enumerate(P):
    r = recs[i]; start = r['parentObj'] if kind == 'E155 parent' else r['childObj']
    for side, st, ms in (('child', start, mv(cm)), ('lna%d' % k, classes[base][k], mv(lm))):
        for lab, a, v, c in edges_of(st, ms):
            TE = TiltEnd(a, v); A, d, carrows, crels = child_info(c)
            for j in range(len(crels)):
                new = crels[:j] + crels[j + 1:]; qd = quotient_dims(carrows, new); extra = sum(qd.values()) - sum(d.values())
                T['edge-drop', 'dim+%d' % extra if extra else 'redundant', 'rejected' if dimcheck(carrows, new, d) else 'accepted'] += 1
for k, x in sorted(T.items()): print(k, x)
