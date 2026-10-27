"""Round 090 (toolsmith), T10/T5: find the children dropped for a `None` canonicalKey in the 7 + 7 J = 0 join of the F-010 pair (E-219/E-226 line, join97).
The join97 README says 26 353 / 15 401 dropped; the log prints 'expansions N None-key children dropped M', so those are the EXPANSION counts and M = 9 / 11.
  toolsmith_nonekey.py scan CK SHARD K BUDGET_S OUT   CK = skeptic_f010join.py `side ... 6` checkpoint (completed depth 6, idx 0); expands frontier[SHARD::K]
        (the depth-6 parents, children = depth 7), records every child whose canonicalKey is None; resumable (OUT is a pickle), exit when BUDGET_S spent.
  toolsmith_nonekey.py report OUTA... -- OUTB...     dedupe by isomorphism, bundle shapes, cap cost, A-vs-B isomorphism test.
Run from the repository root; needs the join97 skeptic_f010join.py preamble (rounds 049/050/079)."""
import sys, os, time, pickle
MYARGV_ = list(sys.argv); sys.path.insert(0, '.')
sys.argv = ['x', 'none']
src = open('workshop/rounds/079/skeptic_f010join.py').read().split("if MODE == 'side':")[0]
src = src.replace("A_ = sys.argv; MODE = A_[1]", "A_ = ['x', 'none']; MODE = 'none'")
exec(compile(src, 'skel', 'exec'))
from quivermutation import fingerprint
def bundleShape(a):
    b = fingerprint.arrowBundles(a.quiver); return sorted(len(v) for v in b.values())
if MYARGV_[1] == 'scan':
    ck, sh, K, budget, out = MYARGV_[2], int(MYARGV_[3]), int(MYARGV_[4]), float(MYARGV_[5]), MYARGV_[6]
    st = pickle.load(open(ck, 'rb')); assert st['d'] == 6 and st['idx'] == 0, (st['d'], st['idx'])
    fr = st['frontier']; mine = list(range(sh, len(fr), K))
    res = pickle.load(open(out, 'rb')) if os.path.exists(out) else dict(pos=0, found=[], nchild=0, secs=0.0, total=len(mine))
    T = time.time()
    while res['pos'] < len(mine) and time.time() - T < budget:
        j = mine[res['pos']]; x = st['seen'][fr[j]][0]
        for kind, v, c in moves(x):
            res['nchild'] += 1
            if fingerprint.canonicalKey(c) is None: res['found'].append((j, kind, v, c))
        res['pos'] += 1
    res['secs'] += time.time() - T
    pickle.dump(res, open(out + '.tmp', 'wb')); os.replace(out + '.tmp', out)
    print('shard', sh, 'of', K, 'parents', res['pos'], '/', len(mine), 'children', res['nchild'], 'None-key', len(res['found']), 'secs %.0f' % res['secs'], flush=True)
elif MYARGV_[1] == 'report':
    i = MYARGV_.index('--'); FA, FB = MYARGV_[2:i], MYARGV_[i + 1:]
    def load(fs):
        L = []
        for f in fs:
            r = pickle.load(open(f, 'rb')); print(f, 'parents', r['pos'], '/', r['total'], 'children', r['nchild'], 'None', len(r['found'])); L += r['found']
        return L
    LA, LB = load(FA), load(FB)
    for tag, L in (('A', LA), ('B', LB)):
        classes = []
        for (j, kind, v, c) in L:
            for cl in classes:
                if sig(cl[0][3]) == sig(c) and iso(cl[0][3], c): cl.append((j, kind, v, c)); break
            else: classes.append([(j, kind, v, c)])
        print(tag, 'None-key children', len(L), 'iso classes', len(classes))
        for cl in classes:
            c = cl[0][3]; print('  class size', len(cl), 'bundles', bundleShape(c), 'arrows', c.quiver.number_of_edges(), 'dim-free; parents', [m[0] for m in cl], 'moves', [(m[1], m[2]) for m in cl])
        globals()['C' + tag] = classes
    hits = [(i, j) for i, a in enumerate(CA) for j, b in enumerate(CB) if sig(a[0][3]) == sig(b[0][3]) and iso(a[0][3], b[0][3])]
    print('A-class vs B-class isomorphic pairs (meets at total 14 the run could not see):', hits)
