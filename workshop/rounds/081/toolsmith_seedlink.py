"""Round 081 (toolsmith), T10: link seed 1 and seed 8 of classes[base] (n = 7, class 2; E-222 rep 2: child ends at seed 1, parent at seed 8)
by J = 0 tiltingPlus F/R chains (E-212 moves), Hom-verified (E-212 horn) with sign-aware iso at the meet.
  toolsmith_seedlink.py REP DA DB     (repo root; REP=2; balls of depth DA from seed 1 and DB from seed 8, separately; prints first meet)"""
import sys, time
REP, DA, DB = (int(x) for x in sys.argv[1:4])
CLS_OF = {0: 1, 1: 2, 2: 2, 3: 2}
sys.argv = ['x', 'none', str(CLS_OF[REP]), '0']
cert = open('workshop/rounds/073/toolsmith_certificate.py').read().split("recs = [r for r")[0]
exec(compile(cert, 'cert', 'exec'))
certIso = iso
sys.argv = ['x', str(REP), '4', '5', '400000']
j = open('workshop/rounds/075/toolsmith_join.py').read().split("# sanity:")[0]
exec(compile(j, 'join', 'exec'))
ij = open('workshop/rounds/075/toolsmith_isojoin.py').read()
ij = ij[ij.index("from quivermutation import arrowPaths"):ij.index("def ballH")]
exec(compile(ij, 'ij', 'exec'))
signedIso = iso; iso = certIso
def run_side(x0, mv):
    ok = True; det = []; x = x0
    for s_ in mv:
        kind, v = s_[0], int(s_[1:]); h = horn(x, kind, v)
        if h is None: return False, ['NOSTEP'], None
        y, zero, cart, _ = h
        det.append('%s%d%s' % (kind, v, '' if zero and cart else '!'))
        if not (zero and cart): ok = False
        x = y
    return ok, det, x
def pball(seed, D):
    k0 = fingerprint.canonicalKey(seed); tab = {k0: (None, None, None)}; objs = {k0: seed}; lvl = {k0: 0}; fr = [(k0, seed)]; hs = {}
    hs.setdefault(H(seed), []).append(k0)
    for d in range(1, D + 1):
        nx_ = []
        for pk, a in fr:
            for kind, v, c in moves(a):
                k = fingerprint.canonicalKey(c)
                if k is None or k in tab: continue
                tab[k] = (pk, kind, v); objs[k] = c; lvl[k] = d; nx_.append((k, c)); hs.setdefault(H(c), []).append(k)
        fr = nx_; print('  depth', d, 'new', len(nx_), flush=True)
    return tab, objs, lvl, hs
def path(tab, k):
    mv = []
    while tab[k][0] is not None: pk, kind, v = tab[k]; mv.append('%s%d' % (kind, v)); k = pk
    return mv[::-1]
S = classes[base]; print('seeds', len(S), 'rep', REP)
kidx = {fingerprint.canonicalKey(s): i for i, s in enumerate(S)}
print('seed keys', [(i, fingerprint.canonicalKey(s) is not None) for i, s in enumerate(S)][:3], 'seed1 key==seed8 key', fingerprint.canonicalKey(S[1]) == fingerprint.canonicalKey(S[8]))
T0 = time.time()
A = pball(S[1], DA); print('ball seed1', len(A[0]), '%.0fs' % (time.time() - T0), flush=True)
B = pball(S[8], DB); print('ball seed8', len(B[0]), '%.0fs' % (time.time() - T0), flush=True)
best = None
for h, ks in A[3].items():
    for ka in ks:
        for kb in B[3].get(h, []):
            tot = A[2][ka] + B[2][kb]
            if (best is None or tot < best[0]) and signedIso(A[1][ka], B[1][kb]): best = (tot, ka, kb)
if not best: print('NO MEET at', DA, '+', DB)
else:
    tot, ka, kb = best; mA, mB = path(A[0], ka), path(B[0], kb)
    okA, dA, xA = run_side(S[1], mA); okB, dB, xB = run_side(S[8], mB)
    meet = signedIso(xA, xB) if xA is not None and xB is not None else False
    print('MEET total', tot, 'seed1 moves', ' '.join(mA) or '-', '| seed8 moves', ' '.join(mB) or '-', '| hom', ' '.join(dA), '||', ' '.join(dB), '| edges-pass', okA and okB, 'signedIso', meet, 'E-212 iso', certIso(xA, xB) if xA is not None and xB is not None else False, 'VERIFIED' if okA and okB and meet else 'NOT-VERIFIED')
