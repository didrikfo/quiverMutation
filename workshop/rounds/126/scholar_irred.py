"""Round 126 (scholar), H-027: is every step of every shortest E-287 tilting-exchange path B_n -> A_n an IRREDUCIBLE silting mutation
(Aihara-Iyama 1009.3370 v3, Def 2.34), not a composite of several mutations at the same summand?

Theory used (AI v3, numbers by counter): for T = D + U basic silting in a Krull-Schmidt T with thick(D) functorially finite, the silting
objects containing D form one orbit under iterated irreducible mutation at the free summand (Thm 2.44, via silting reduction Thm 2.37,
compatibility Lemma 2.40 and Thm 2.26 silt U = {M[i]}); so they are mu^k(T), k in Z, a chain T > mu^+T > mu^+mu^+T > ... (Thm 2.35).
Hence for T' = D + y in that fibre, T' is an irreducible mutation of T iff T, T' are comparable and NO silting D + Z lies strictly between.
Test: count indecomposables Z (all n(n-1) objects of D^b(B_n) ~ D^b(kD_n) up to shift, shifts in [-W, W]), Z not U, y, not in D, with
D + Z presilting (Hom(.,.[i]) = 0, i > 0), K0-unimodular, and T >= D + Z >= T' (left) or T' >= D + Z >= T (right).  Every silting object
is presilting and unimodular (AI Thm 2.27), so count 0 => irreducible (sound); count > 0 => not certified (the candidate could fail to
generate), reported with the candidates.
Controls: (1) all depth-1 exchanges of B_n: at most 2 complements per slot can be irreducible (mu^+ and mu^-), so a slot with 3 tilting
complements must give count > 0 at least once; (2) planted composite: for each pair of consecutive depth-1/2 exchanges at the SAME slot in
the same direction (T > T1 > T2, all tilting), the direct exchange T -> T2 must give count >= 1.
Usage: scholar_irred.py n    (repository root; n = 5, 6, 7)"""
import sys, itertools, json, os, time
sys.path.insert(0, '.')
exec(open('workshop/rounds/107/toolsmith_tiltsearch.py').read().split("\nif __name__ == '__main__':")[0])

n = int(sys.argv[1]); W = int(os.environ.get('W', 7)); t0 = time.time()   # shift window; W=10 re-run as a window check
alg = Bn(n); A = Alg(alg)
f = 'workshop/rounds/107/toolsmith_objs_%d.json' % n
if not os.path.exists(f): f = 'workshop/rounds/109/toolsmith_objs_%d.json' % n
objs = [deser(S) for S in json.load(open(f))]; N = len(objs); J = 14
h = [[{j: homdim(A, objs[r], shiftC(objs[s], j), 0) for j in range(-J, J + 1)} for s in range(N)] for r in range(N)]
print('n', n, 'objects up to shift', N, '(expected n(n-1) =', n * (n - 1), ') Hom tables %.0fs' % (time.time() - t0)); sys.stdout.flush()
cls = [np.array(klassC(X, n)) for X in objs]
def hom(x, y, m=0):
    (r, a), (s, b) = x, y; j = b - a + m
    return h[r][s][j] if -J <= j <= J else 0
POS = range(1, 2 * W + 3); NONZ = [m for m in range(-2 * W - 2, 2 * W + 3) if m != 0]
cache = {}
def cp(x, y):
    if (x, y) not in cache: cache[(x, y)] = all(hom(x, y, m) == 0 and hom(y, x, m) == 0 for m in NONZ)
    return cache[(x, y)]
selfok = lambda x: all(hom(x, x, m) == 0 for m in range(-8, 9) if m != 0) and hom(x, x, 0) == 1
klass = lambda x: ((-1) ** x[1]) * cls[x[0]]
def unimod(T): return abs(round(np.linalg.det(np.array([klass(x) for x in T])))) == 1
def geq(S, T):   # S >= T in AI's order: Hom(S, T[i]) = 0 for all i > 0
    return all(hom(x, y, m) == 0 for x in S for y in T for m in POS)
def presilt(S): return geq(S, S)
P = {}
for r, X in enumerate(objs):
    if len(X[0]) == 1 and len(X[0][0]) == 1: P[X[0][0][0]] = (r, 0)
T0 = tuple(P[i] for i in range(1, n + 1))
allz = [(r, a) for r in range(N) for a in range(-W, W + 1)]
CA = np.array(invariants.cartanMatrix(An(n), exact=True).tolist(), dtype=int)
def Hm(T): return np.array([[hom(x, y) for y in T] for x in T])
def isA(T):
    M = Hm(T)
    return sorted(M.flatten()) == sorted(CA.flatten()) and (perm_canon_eq(M, CA) or perm_canon_eq(M.T, CA))
def nbrs(T):
    out = []
    for i, U in enumerate(T):
        rest = [x for x in T if x != U]
        for y in allz:
            if y == U or y in rest or not selfok(y): continue
            if all(cp(y, x) for x in rest):
                T2 = T[:i] + (y,) + T[i + 1:]
                if unimod(T2): out.append((i, T2))
    return out

def irred(T, T2):
    """returns (direction, count of presilting unimodular D+Z strictly between, the Zs)."""
    D = [x for x in T if x in T2]; assert len(D) == n - 1
    U = [x for x in T if x not in T2][0]; y = [x for x in T2 if x not in T][0]
    left, right = geq(T, T2), geq(T2, T)
    if left == right: return ('incomparable' if not left else 'both'), None, []
    hi, lo = (T, T2) if left else (T2, T)
    Zs = []
    for z in allz:
        if z in (U, y) or z in D or not selfok(z): continue
        S = D + [z]
        if presilt(S) and unimod(S) and geq(hi, S) and geq(S, lo): Zs.append(z)
    return ('left' if left else 'right'), len(Zs), Zs

# ---- control 1 and 2: depth-1 and depth-2 exchanges of B_n
print('== control 1: depth-1 exchanges of B_%d (slot, new summand, direction, #between) ==' % n)
d1 = nbrs(T0); byslot = {}
for i, T1 in d1:
    dr, c, Zs = irred(T0, T1); byslot.setdefault(i, []).append(c)
    print('  slot %d: %s -> %s  %s  between=%s %s' % (i, T0[i], T1[i], dr, c, Zs))
bad = [i for i, cs in byslot.items() if sum(1 for c in cs if c == 0) > 2]
print('  slots with > 2 certified-irreducible complements (impossible, would be a test failure):', bad)
print('  slots with 3+ complements, of which the test says NOT irreducible:', {i: sum(1 for c in cs if c) for i, cs in byslot.items() if len(cs) >= 3})
print('== control 2: planted composites T0 > T1 > T2 (or <), same slot, same direction, all tilting; direct T0 -> T2 must count >= 1 ==')
planted = caught = 0
for i, T1 in d1:
    dr1, c1, _ = irred(T0, T1)
    if c1 != 0: continue
    for j, T2 in nbrs(T1):
        if j != i or T2[i] == T0[i]: continue
        dr2, c2, _ = irred(T1, T2)
        if c2 != 0 or dr2 != dr1: continue
        dr, c, Zs = irred(T0, T2); planted += 1; caught += (c is not None and c >= 1)
        print('  slot %d: %s -> %s -> %s (%s,%s): direct %s between=%s %s' % (i, T0[i], T1[i], T2[i], dr1, dr2, dr, c, Zs))
print('  planted composites caught: %d of %d' % (caught, planted))
sys.stdout.flush()

# ---- all shortest paths to the first-depth A_n hits (as round 109)
key_ = frozenset
layers = [{key_(T0): T0}]; seen = {key_(T0)}; parents = {key_(T0): []}
for dpt in range(1, 7):
    new = {}
    for ks, T in layers[-1].items():
        for i, T2 in nbrs(T):
            k2 = key_(T2)
            if k2 in seen and k2 not in new: continue
            new.setdefault(k2, T2); parents.setdefault(k2, []).append((ks, T, T2))
    seen |= set(new); layers.append(new)
    hits = [k for k in new if isA(new[k])]
    print('depth', dpt, 'new', len(new), 'A_n hits', len(hits), '[%.0fs]' % (time.time() - t0)); sys.stdout.flush()
    if hits: break
order = {}
for d, L in enumerate(layers):
    for k in L: order[k] = L[k]
def paths(k):
    if not parents[k]: return [[]]
    out = []
    for pk, T, T2 in parents[k]:
        for p in paths(pk): out.append(p + [(T, T2)])
    return out
steps = {}; npaths = 0
for k in hits:
    for p in paths(k):
        npaths += 1
        for T, T2 in p: steps[(key_(T), key_(T2))] = (T, T2)
print('== E-287 shortest paths: %d paths, %d distinct steps ==' % (npaths, len(steps)))
tally = {}
for (a, b), (T, T2) in sorted(steps.items(), key=lambda kv: (sorted(kv[0][0]), sorted(kv[0][1]))):
    dr, c, Zs = irred(T, T2)
    i = [s for s in range(n) if T[s] != T2[s]][0]
    tally[(dr, c == 0)] = tally.get((dr, c == 0), 0) + 1
    print('  %s -> %s at slot %d: %s, between=%s %s' % (T[i], T2[i], i, dr, c, Zs))
print('tally (direction, certified irreducible): count', tally)
# per path: all steps certified?
allok = 0
for k in hits:
    for p in paths(k):
        allok += all(irred(T, T2)[1] == 0 for T, T2 in p)
print('paths all of whose steps are certified irreducible: %d of %d   [%.0fs]' % (allok, npaths, time.time() - t0))
mult = {}
for k in hits:
    for p in paths(k):
        for T, T2 in p:
            dr = irred(T, T2)[0]; mult[dr] = mult.get(dr, 0) + 1
print('path-steps with multiplicity by AI direction (left = mu^+, T > T\'):', mult, '(compare E-289 gate-matching side: right 9/16/25, left 3/4/5)')
