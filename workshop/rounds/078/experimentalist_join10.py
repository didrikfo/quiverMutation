"""S-1 / round 078: tilting-only (J = 0, E-212/E-049 F and R moves, key kept) join search for the equivalent pair
90000000 / 00000030 at n = 10 (rows have n - 2 = 8 entries).  Moves, canonicalKey = rounds/049/toolsmith_tiltpath.py (exec'd, base key reset).
  experimentalist_join10.py plan                 key-class sizes, branching, per-level timing (small depths)
  experimentalist_join10.py run DA DB CAP [F]    alternate BFS from both ends (levels DA, DB), stop at first shared canonicalKey; cap = max keys per side
Run from repo root.  A MISS is a bound for these depths/caps, not a proof of inequivalence (E-164 relabelling effect)."""
import sys, time
MYARGV = sys.argv; sys.argv = ['x', 'none', '0', '1', '1', '1', 'x']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
from quivermutation import coxeterTables as ct
N = 10
A_ROW, B_ROW = (9,0,0,0,0,0,0,0), (0,0,0,0,0,0,3,0)
def alg(r): return nk.LinearNakayamaAlgebra(N, list(r))
A, B = alg(A_ROW), alg(B_ROW)
kA, kB = search._coxeterKeyOrNone(A), search._coxeterKeyOrNone(B)
base = kA
print('keys equal', kA == kB, 'key', str(kA)[:80], flush=True)
def sideBall(start, depth, cap, label):
    k0 = fingerprint.canonicalKey(start); seen = {k0: 0}; fr = [start]; lv = [1]; t = time.time()
    for d in range(1, depth + 1):
        nxt = []
        for a in fr:
            if len(seen) >= cap: break
            for kind, v, c in moves(a):
                k = fingerprint.canonicalKey(c)
                if k in seen: continue
                seen[k] = d; nxt.append(c)
        lv.append(len(nxt)); fr = nxt
        print('  %s level %d: %d new, %d keys, %.0fs' % (label, d, len(nxt), len(seen), time.time() - t), flush=True)
        if len(seen) >= cap or not nxt: break
    return seen, lv
if MYARGV[1] == 'plan':
    t = time.time()
    rows = list(nk.allRelationLengths(N)); same = 0
    for r in rows:
        if ct.lnaCoxeterKey(N, tuple(r)) == kA: same += 1
    print('n=10 rows', len(rows), 'rows with key of 90000000:', same, '%.0fs' % (time.time() - t))
    for nm, a in (('A=90000000', A), ('B=00000030', B)):
        t = time.time(); ms = moves(a); print(nm, 'moves at root', [(k, v) for k, v, c in ms], '%.1fs' % (time.time() - t), flush=True)
        sideBall(a, 4, 3000, nm)
    print('stats', stats)
else:
    DA, DB, CAP = int(MYARGV[2]), int(MYARGV[3]), int(MYARGV[4])
    sa, la = sideBall(A, DA, CAP, 'A'); sb, lb = sideBall(B, DB, CAP, 'B')
    common = set(sa) & set(sb)
    print('A levels', la, 'B levels', lb, 'shared keys', len(common))
    if common: print('HIT total depth', min(sa[k] + sb[k] for k in common))
    else: print('MISS (bound for depth', DA, '+', DB, 'cap', CAP, ')')
