"""Round 078 (toolsmith), T10: record the F/R move lists of the 8 E-216 chains (4 children, 4 parents; n = 7 classes 1, 2) and Hom-verify them.
  toolsmith_chains.py REP DT DC      (repo root; REP = 0..3 as rounds/075/toolsmith_join.py; DT = 4, DC = 5 as 075)
Search = rounds/075/toolsmith_isojoin.py (same moves, WL hash + sign-aware iso at meets) but with parent pointers so the move lists are kept.
Verification of each chain  X --moves--> M  ~iso  M' <--moves-- seed:  every edge by the E-212 Hom test (rounds/073 `horn`: Hom(T,T[-1]) = Hom(T,T[1]) = 0 and
dim Hom(T_i,T_j) = Cartan of the child; no gate/tiltingPlus/key), the meet by the sign-aware iso AND the (necessary) E-212 `iso`.  Generation assumed (E-168)."""
import sys, time
REP, DT, DC = (int(x) for x in sys.argv[1:4])
CLS_OF = {0: 1, 1: 2, 2: 2, 3: 2}
ns = globals()
sys.argv = ['x', 'none', str(CLS_OF[REP]), '0']
cert = open('workshop/rounds/073/toolsmith_certificate.py').read().split("recs = [r for r")[0]
exec(compile(cert, 'cert', 'exec'))
certIso = iso          # E-212 necessary test
certStruct = struct
sys.argv = ['x', str(REP), str(DT), str(DC), '400000']
j = open('workshop/rounds/075/toolsmith_join.py').read().split("# sanity:")[0]
exec(compile(j, 'join', 'exec'))
ij = open('workshop/rounds/075/toolsmith_isojoin.py').read()
ij = ij[ij.index("from quivermutation import arrowPaths"):ij.index("def ballH")]
exec(compile(ij, 'ij', 'exec'))     # struct, H, iso (sign-aware)
signedIso = iso; iso = certIso
def run_side(x0, mv, _):
    ok = True; ne = 0; det = []; x = x0
    for s_ in mv:
        kind, v = s_[0], int(s_[1:]); h = horn(x, kind, v)
        if h is None: return False, ne, ['NOSTEP'], None
        y, zero, cart, hm1 = h; ne += 1
        det.append('%s%d%s' % (kind, v, '' if zero and cart else '!'))
        if not (zero and cart): ok = False
        x = y
    return ok, ne, det, x
T0 = time.time()
def pathTo(tab, k):
    mv = []
    while tab[k][1] is not None:
        _, pk, kind, v = tab[k]; mv.append('%s%d' % (kind, v)); k = pk
    return mv[::-1]
# target ball with parents
tab = {}; tobjs = {}; frontier = []
for s in classes[base]:
    k = fingerprint.canonicalKey(s)
    if k in tab: continue
    tab[k] = (0, None, None, None); tobjs.setdefault(H(s), []).append((s, k, 0)); frontier.append((k, s))
for d in range(1, DT + 1):
    nxt = []
    for pk, alg in frontier:
        for kind, v, c in moves(alg):
            k = fingerprint.canonicalKey(c)
            if k is None or k in tab: continue
            tab[k] = (d, pk, kind, v); tobjs.setdefault(H(c), []).append((c, k, d)); nxt.append((k, c))
    frontier = nxt
print('rep', REP, 'class', CLS, 'target ball keys', len(tab), '%.0fs' % (time.time() - T0), flush=True)
def seedOf(k):
    while tab[k][1] is not None: k = tab[k][1]
    return k
seedObj = {fingerprint.canonicalKey(s): s for s in classes[base]}
seedIdx = {}
for _i, _s in enumerate(classes[base]): seedIdx.setdefault(fingerprint.canonicalKey(_s), _i)   # r078 review item 3: index of the seed in classes[base]
def findMeet(alg):
    k0 = fingerprint.canonicalKey(alg); seen = {k0: (None, None, None)}; objs = {k0: alg}; frontier = [(k0, alg)]; best = None
    def check(k, a, d):
        nonlocal best
        for (o, ok, td) in tobjs.get(H(a), []):
            if (best is None or d + td < best[0]) and signedIso(a, o): best = (d + td, k, ok, d, td)
    check(k0, alg, 0)
    for d in range(1, DC + 1):
        nxt = []
        for pk, a in frontier:
            for kind, v, c in moves(a):
                k = fingerprint.canonicalKey(c)
                if k is None or k in seen: continue
                seen[k] = (pk, kind, v); objs[k] = c; nxt.append((k, c)); check(k, c, d)
        frontier = nxt
        if best: break
    tot, ak, tk, d, td = best
    mv = []; k = ak
    while seen[k][0] is not None: pk, kind, v = seen[k]; mv.append('%s%d' % (kind, v)); k = pk
    return mv[::-1], objs[ak], pathTo(tab, tk), seedObj[seedOf(tk)], seedObj.get(seedOf(tk))
for name, alg in (('child', C), ('parent', P)):
    mA, meetA, mB, seed, _ = findMeet(alg)
    okA, nA, dA, xA = run_side(alg, mA, None)
    okB, nB, dB, xB = run_side(seed, mB, None)
    meet = signedIso(xA, xB) if xA is not None and xB is not None else False
    meet2 = certIso(xA, xB) if xA is not None and xB is not None else False
    print(name, 'rep', REP, 'moves(alg->meet)', ' '.join(mA) or '-', '| moves(seed->meet)', ' '.join(mB) or '-', '| seed key', fingerprint.canonicalKey(seed) in tab and 'in ball', '| seed index', seedIdx[fingerprint.canonicalKey(seed)], 'of', len(classes[base]))
    print(name, 'rep', REP, 'edges', nA + nB, 'hom-detail', ' '.join(dA), '||', ' '.join(dB), '| all-edges-pass', okA and okB, '| meet signedIso', meet, 'E-212 iso', meet2,
          'VERIFIED' if (okA and okB and meet) else 'NOT-VERIFIED', flush=True)
