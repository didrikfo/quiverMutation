"""Round 057 (toolsmith), T10 item 1.  Usage (repo root; needs /tmp/tsm/c<CLS>.pkl from
  DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 CLS 20000 /tmp/tsm/cCLS.pkl 100   (rerun until it finishes)):
  CLS=1 toolsmith_endt2_run.py fail                 End(T) vs child at every failing J != 0 key-keeping step of class CLS (matrix-valued arrows)
  CLS=1 toolsmith_endt2_run.py plan                 list the E-155 / E-158 paths of class CLS with edge counts (no End(T))
  CLS=1 toolsmith_endt2_run.py paths LO:HI          End(T) on both sides of paths LO..HI-1 of the plan list (E-155 child, E-155 parent, E-158)
"""
import sys, re, time, pickle
MODE = sys.argv[1]; SL = sys.argv[2] if len(sys.argv) > 2 else None
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
E158 = {2: [('E158 child', 6, 'F1 F7 R3 R5 R5 F4 R6', 0, 'F1 F1 F2 R6 R5')],
        1: [('E158 child', 5, 'F1 F3 F5 F4 R7 F3 R1', 9, 'F7 R2 R1 F2 F7'), ('E158 child', 12, 'R7 R4 R7 R1', 6, 'F1 F5 F4 R7 F3'),
            ('E161 child', 14, 'F1 F3 F1 F5 F4 R7 R1', 9, 'F7 R2 R1 F2 F7 R3')]}

def plan():
    out = []; txt = open('workshop/rounds/049/toolsmith_paths_logs.txt').read().split('== ')
    for sec in txt[1:]:
        name = sec.split('\n')[0].strip()
        if not name.endswith('_c%d' % CLS): continue
        kind = 'E155 child' if name.startswith('paths') else 'E155 parent'
        for ln in sec.split('\n'):
            m = re.match(r'(paths|parents) (\d+) key \w+ parentdepth \d+ replay ok total \d+ \| child->meet: (.*?) \| LNA/dual #(\d+) ->meet: (.*?) \d+s', ln)
            if m: out.append((kind, int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)))
    return out + E158[CLS]

def mv(s): return [] if s.strip() == '-' else s.split()

if MODE == 'fail':
    for n_, r in enumerate(recs):
        a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v)); t0 = time.time()
        rr = compare2(a, v, c); print('fail c%d' % CLS, n_, 'v', v, 'dims', rr['dims'], 'arrows', rr['arrows'], rr['n_arrows'], 'maxpar', rr['par'], '|', rr['verdict'], '%.0fs' % (time.time() - t0), flush=True)
    sys.exit()
P = plan()
if MODE == 'plan':
    for n_, (kind, i, cm, k, lm) in enumerate(P): print(n_, kind, 'idx', i, 'edges', len(mv(cm)) + len(mv(lm)), '|', cm, '|', 'LNA#%d' % k, lm)
    print('paths', len(P), 'edges', sum(len(mv(p[2])) + len(mv(p[4])) for p in P)); sys.exit()
lo, hi = (int(x) for x in SL.split(':')); tot = {}
for n_, (kind, i, cm, k, lm) in enumerate(P):
    if not lo <= n_ < hi: continue
    r = recs[i]; start = r['parentObj'] if kind == 'E155 parent' else r['childObj']
    for side, st, ms in (('child', start, mv(cm)), ('lna%d' % k, classes[base][k], mv(lm))):
        for lab, a, v, c in edges_of(st, ms):
            t0 = time.time(); rr = compare2(a, v, c); vd = rr['verdict'].split(' (')[0]
            tot[vd] = tot.get(vd, 0) + 1
            print('path', n_, kind, i, side, lab, 'dims', rr['dims'], 'arrows', rr['arrows'], 'par', rr['par'], '|', rr['verdict'], '%.1fs' % (time.time() - t0), flush=True)
print('TOTAL', tot)
