"""Round 081 (toolsmith), response to referee items 2-4: replay the E-222 rep 2 chain child -> seed 1 -> seed 8 -> parent (2 + 4 + 5 = 11 edges)
as ONE run in ONE process on the same objects, with per-edge Hom dimension lines, plus a single-move corruption control on the 4 seed-link edges.
  toolsmith_chain11.py      (repo root; ~3-4 min)
Edges are verified in the direction seed -> end (E-212 `horn` is directional): P1 seed1 -F6 F7-> xc ~ child; P2 seed1 -F6-> M ~ M' <-R5 R5 R6- seed8;
P3 seed8 -F1 F1 F2 R7-> xp ~ yp <-F5- parent.  Meets by sign-aware iso AND E-212 `iso`.  The three pieces share the same seed objects (S[1], S[8])."""
import sys, time
REP = 2; CLS_OF = {0: 1, 1: 2, 2: 2, 3: 2}
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
S = classes[base]; T0 = time.time()
def edge(x, kind, v):
    r = step(x, kind, v)
    if r is None: return None
    y, a, c = r
    Vs, res = stepTest(a, v, -1)
    m1 = np.array(res[-1]); p1 = np.array(res[1]); H0 = np.array(res[0], dtype=int)
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    zero = not m1.any() and not p1.any(); cart = bool((H0 == C.T).all())
    return y, zero, cart, int(np.abs(m1).sum()), int(np.abs(p1).sum()), int(H0.sum()), int(C.sum())
def run(x0, mv, label, quiet=False):
    ok = True; x = x0; lines = []
    for i, s_ in enumerate(mv):
        kind, v = s_[0], int(s_[1:]); h = edge(x, kind, v)
        if h is None:
            lines.append('  %s edge %d %s NOSTEP' % (label, i + 1, s_)); return False, lines, None
        y, zero, cart, m1, p1, h0, c0 = h
        lines.append('  %s edge %d %s%d: sum Hom(T,T[-1]) = %d, sum Hom(T,T[+1]) = %d, sum dim Hom(T,T) = %d, sum Cartan(child) = %d, Cartan match %s -> %s' % (label, i + 1, kind, v, m1, p1, h0, c0, cart, 'pass' if zero and cart else 'REJECT'))
        if not (zero and cart): ok = False
        x = y
    return ok, lines, x
def both(x, y): return signedIso(x, y), certIso(x, y)
allok = True; nedge = 0
def piece(name, pairs, meet_with):
    global allok, nedge
    print(name)
    ends = []
    for label, x0, mv in pairs:
        ok, lines, x = run(x0, mv, label); print('\n'.join(lines)); nedge += len(mv); allok &= ok; ends.append(x)
    return ends
# P1: seed 1 -F6 F7-> xc ~ child
(xc,) = piece('P1 child <-> seed 1', [('seed1', S[1], ['F6', 'F7'])], None)
m1 = both(xc, C); print('  meet P1: end of seed1 F6 F7 vs CHILD: signedIso %s, E-212 iso %s' % m1); allok &= m1[0]
# P2: seed link
xa, xb = piece('P2 seed 1 <-> seed 8', [('seed1', S[1], ['F6']), ('seed8', S[8], ['R5', 'R5', 'R6'])], None)
m2 = both(xa, xb); print('  meet P2: signedIso %s, E-212 iso %s' % m2); allok &= m2[0]
print('  seed1 != seed8: signedIso(S1,S8) = %s' % signedIso(S[1], S[8]))
# P3: seed 8 -F1 F1 F2 R7-> xp ~ yp <-F5- parent
xp, yp = piece('P3 seed 8 <-> parent', [('seed8', S[8], ['F1', 'F1', 'F2', 'R7']), ('parent', P, ['F5'])], None)
m3 = both(xp, yp); print('  meet P3: signedIso %s, E-212 iso %s' % m3); allok &= m3[0]
print('CHAIN child ~ xc <-F7 F6- seed1 -F6-> M ~ M\' <-R5 R5 R6- seed8 -F1 F1 F2 R7-> xp ~ yp <-F5- parent: edges %d, all Hom+Cartan pass and 3 meets signed-iso: %s  (%.0fs)' % (nedge, allok, time.time() - T0), 'VERIFIED' if allok and nedge == 11 else 'NOT-VERIFIED', flush=True)
# corruption control on the 4 seed-link edges: replace one move by every other (kind, v)
print('CORRUPTION CONTROL (seed link: seed1 F6 || seed8 R5 R5 R6; one move replaced by each other (kind, v), v in 1..7)')
V = sorted(S[1].quiver.nodes); tot = nostep = hrej = pass_edges = pass_all = 0
for side, x0, mv in (('seed1', S[1], ['F6']), ('seed8', S[8], ['R5', 'R5', 'R6'])):
    for i in range(len(mv)):
        for kind in 'FR':
            for v in V:
                if '%s%d' % (kind, v) == mv[i]: continue
                mv2 = list(mv); mv2[i] = '%s%d' % (kind, v); tot += 1
                ok, lines, x = run(x0, mv2, side)
                if x is None and not ok and 'NOSTEP' in lines[-1]: nostep += 1; continue
                if not ok: hrej += 1; continue
                pass_edges += 1
                # edges pass; does the meet survive?
                other = (xb if side == 'seed1' else xa)
                if signedIso(x, other): pass_all += 1
print('  corruptions %d: NOSTEP %d, Hom/Cartan REJECT %d, all edges pass %d, of those meet still signed-iso %d  (%.0fs)' % (tot, nostep, hrej, pass_edges, pass_all, time.time() - T0))
