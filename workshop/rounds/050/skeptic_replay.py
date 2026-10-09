"""Round 050 (skeptic), T10 (i): replay the printed E-157 paths (rounds/049/toolsmith_paths_logs.txt) and test EVERY edge with the independent
tilting test of skeptic_tilt.py (no tiltingPlus / perI / gate used for the verdict; they are printed for context).
Usage (repo root): skeptic_replay.py in.pkl cls   where in.pkl = rounds/049/toolsmith_collect.py output for class cls (7 1 / 7 2).
Per edge x -> y (move F v: step at v of x; move R v: step at v of dual(x), carried back): a = x or dual(x), c = reduced mutation of a at v.
Verdict 'tilt' := Hom(T,T[-1]) = 0 and Hom(T,T[1]) = 0 for T = T_v(a) of skeptic_tilt.py; 'cartan' := H = Cartan(c)^T (labelled, c = child before dualizing).
Also checks the paths: child-side end key == LNA-side end key."""
import sys, re, pickle
import numpy as np
A_ = sys.argv; PK, CLS = A_[1], int(A_[2])
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import hashlib
from quivermutation import invariants
def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]
def step(x, kind, v):
    """returns (y, a, c_a, V) or None"""
    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    V = sorted(x.quiver.nodes)
    raw = mutation.quiverMutationAtVertex(a, v)
    c = reduction.reducePathAlgebra(raw)
    if sorted(c.vertices()) != V: return None
    y = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return y, a, c, V
def verdict(x, kind, v):
    r = step(x, kind, v)
    if r is None: return None
    y, a, c, V = r
    Vs, res = stepTest(a, v, -1)
    zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    H = np.array(res[0], dtype=int)
    cart = bool((H == C.T).all())
    ctx = dict(J=bool(perI(a, v)), tp=bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v)))
    return y, zero, cart, ctx, int(np.abs(np.array(res[-1])).sum())
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
seedkeys = classes[base]
tot = dict(edges=0, tilt=0, cart=0, both=0, J=0, tp=0, paths=0, pathsOK=0, bad=[])
SEC = None
for line in open('workshop/rounds/049/toolsmith_paths_logs.txt'):
    if line.startswith('=='): SEC = line.strip('= \n'); continue
    if SEC is None or not SEC.endswith('_c%d' % CLS): continue
    m = re.match(r'(paths|parents) (\d+) key (\w+) parentdepth (\d+) replay ok total (\d+) \| child->meet: (.*?) \| LNA/dual #(\d+) ->meet: (.*?) \d+s', line)
    if not m: continue
    which, i, key, pd, total, cs, lna, ts = m.groups(); i = int(i); lna = int(lna)
    if len(A_) > 3 and which != A_[3]: continue
    r = recs[i]; x0 = r['childObj'] if which == 'paths' else r['parentObj']
    seqs = [(x0, cs.split()), (seedkeys[lna], ts.split())]
    ends = []; allok = True; out = []
    for x, mv in seqs:
        for s_ in mv:
            kind, v = s_[0], int(s_[1:])
            vd = verdict(x, kind, v)
            if vd is None: allok = False; out.append(s_ + ':NOSTEP'); break
            y, zero, cart, ctx, hm1 = vd
            tot['edges'] += 1; tot['tilt'] += zero; tot['cart'] += cart; tot['both'] += (zero and cart); tot['J'] += ctx['J']; tot['tp'] += ctx['tp']
            if not (zero and cart): tot['bad'].append((which, i, s_, zero, cart, hm1, ctx)); allok = False
            out.append('%s%s%s' % (s_, '' if zero and cart else '!', 'J' if ctx['J'] else ''))
            x = y
        ends.append(fingerprint.canonicalKey(x))
    meet = ends[0] == ends[1]
    tot['paths'] += 1; tot['pathsOK'] += (allok and meet)
    print(which, i, key, 'total', total, 'meet', meet, 'allEdgesTilting', allok, '|', ' '.join(out[:len(cs.split())]), '||', ' '.join(out[len(cs.split()):]), flush=True)
print('SUMMARY cls', CLS, {k: v for k, v in tot.items() if k != 'bad'}); print('BAD', tot['bad'])
