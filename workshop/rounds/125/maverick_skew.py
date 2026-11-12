"""Round 125 maverick, T10 / E-318: a derived invariant that is NOT a function of S = C + C^T alone, on the n = 11 key group of 72
(toolsmith_jxgroups.json['11'][5]: 13 + 59 by SNF(S)), and on every node of E-318's 5-step walk 003033030 -F3 R8 F3 F3 F6-> 030330030.

Invariants (each an invariant of Z-congruence C -> P^T C P, P in GL_n(Z), hence of derived equivalence for finite global dimension
(H-024 necessity step, Ladkani math/0610685 Cor. 3.13 / BHL 0906.3422 Prop. 2.1)):
  SK     SNF(C - C^T)                    (skew part; = coker(Phi + I) up to the unit C^{-1})
  X<a>   SNF(a C + C^T), a in {2, 3, -2, -3}   (pencil; a = 1 is S, a = -1 is SK)
  F<f>   SNF(f(Phi)) for each irreducible factor f of the Coxeter polynomial, Phi = -C^{-1} C^T (conjugation invariant)
Reference (functions of S alone): SNF(S).
Controls:
  census  every move of freeMoves/lnaMoves/edgeMoves/doubleMutation (recorded derived equivalences) inside the group: invariants must be
          constant on each census orbit (YES control on known-equivalent pairs).
  cong    planted Z-congruent copies P^T C P (random unimodular P, 30 per class representative): must be equal (catches non-invariant code).
  walk    the E-318 walk: the J = 0 tilting steps must keep every invariant; the crossing step F3 at P is reported.
SNF by python-flint, cross-checked by pari for every matrix.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skew.py [walk=1]"""
import sys, os, json, random, collections, time
sys.path.insert(0, '.')
import flint, cypari, sympy as sp
from quivermutation import coxeterTables as ct, freeMoves, lnaMoves, edgeMoves, doubleMutation
pari = cypari.pari; N = 11; T0 = time.time()


def snf(M):
    M = [[int(v) for v in row] for row in M]; n = len(M); a = flint.fmpz_mat(M).snf(); d1 = tuple(sorted(abs(int(a[i, i])) for i in range(n)))
    d2 = tuple(sorted(abs(int(x)) for x in pari.matsnf(pari.matrix(n, n, [v for r in M for v in r]))))
    if len(d2) < n: d2 = tuple(sorted((0,) * (n - len(d2)) + d2))
    assert d1 == d2, (d1, d2)
    return tuple(x for x in d1 if x != 1)          # drop the 1s (rank info kept through the 0s)


def invs(C):
    n = len(C); Cm = sp.Matrix(C); T = Cm.T
    out = {'S': snf((Cm + T).tolist()), 'SK': snf((Cm - T).tolist())}
    for a in (2, 3, -2, -3): out['X%d' % a] = snf((a * Cm + T).tolist())
    Phi = -Cm.inv() * T; assert all(x.is_integer for x in Phi)
    x = sp.Symbol('x'); cp = sp.Poly((x * sp.eye(n) - Phi).det(), x)
    for f, _ in sp.factor_list(cp.as_expr())[1]:
        fp = sp.Poly(f, x); M = sp.zeros(n, n); P = sp.eye(n)
        for c in reversed(fp.all_coeffs()): M += c * P; P = P * Phi
        out['F[%s]' % str(f).replace(' ', '')] = snf(M.tolist())
    return out


def lna(rl): return [list(r) for r in ct.lnaCartanMatrix(N, tuple(rl))]


def unimod(n, rng, k=25):
    P = sp.eye(n)
    for _ in range(k):
        i, j = rng.sample(range(n), 2); E = sp.eye(n); E[i, j] = rng.choice([-2, -1, 1, 2]); P = P * E
    return P


def main():
    g = json.load(open('workshop/rounds/122/toolsmith_jxgroups.json'))['11'][5]
    mem = sorted(tuple(int(c) for c in s) for part in g for s in part); memset = set(mem)
    I = {r: invs(lna(r)) for r in mem}; keys = sorted(I[mem[0]])
    print('group of %d, invariants %s  (%.0fs)' % (len(mem), keys, time.time() - T0), flush=True)
    A, B = tuple(int(c) for c in '003033030'), tuple(int(c) for c in '030330030')
    print('\n== the E-318 pair (SNF entries != 1 shown)')
    for k in keys: print('  %-22s 003033030 %-18s 030330030 %-18s %s' % (k, I[A][k], I[B][k], 'SEPARATES' if I[A][k] != I[B][k] else 'equal'))
    # how each invariant partitions the group, and whether it refines / is refined by SNF(S)
    print('\n== partition of the 72 by each invariant: classes, and joint with SNF(S)')
    for k in keys:
        c = collections.Counter(I[r][k] for r in mem); j = collections.Counter((I[r]['S'], I[r][k]) for r in mem)
        print('  %-22s %d classes %s | joint with S: %d classes' % (k, len(c), sorted(c.values(), reverse=True), len(j)))
    full = collections.Counter(tuple(I[r][k] for k in keys) for r in mem)
    print('  all together: %d classes, sizes %s' % (len(full), sorted(full.values(), reverse=True)))
    # census orbits (YES control)
    parent = {r: r for r in mem}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    ne = 0
    for r in mem:
        tg = [tuple(x) for x in lnaMoves.rewritesOf(N, list(r))] + [tuple(x) for x, _ in edgeMoves.rewritesOf(N, r)]
        tg += [tuple(x) for x, _ in doubleMutation.rewritesOf(N, r)]
        s = freeMoves.stripLengthTwo(r); tg += [tuple(s)] if s != r else []
        for v in freeMoves.addableLengthTwo(N, r): d = list(r); d[v - 1] = 2; tg.append(tuple(d))
        for x in tg:
            if x in memset: ne += 1; parent[find(r)] = find(x)
    orb = collections.defaultdict(list)
    for r in mem: orb[find(r)].append(r)
    bad = [(k, o) for o in orb.values() for k in keys if len({I[r][k] for r in o}) > 1]
    print('\n== census control: %d in-group move edges, %d orbits sizes %s; orbits on which some invariant varies: %d (must be 0)'
          % (ne, len(orb), sorted(map(len, orb.values()), reverse=True), len(bad)))
    for o in sorted(orb.values(), key=len, reverse=True):
        print('  orbit %2d  S %-10s SK %-12s  %s' % (len(o), I[o[0]]['S'], I[o[0]]['SK'], ' '.join(ct.className(r) for r in o[:5])))
    # planted congruent copies (YES control)
    rng = random.Random(125); reps = [o[0] for o in orb.values()][:6]; fails = 0; tot = 0
    for r in reps:
        C = sp.Matrix(lna(r))
        for _ in range(30):
            P = unimod(N, rng); assert abs(P.det()) == 1; tot += 1
            if invs((P.T * C * P).tolist()) != I[r]: fails += 1
    print('\n== congruence control: %d planted P^T C P copies, %d with a changed invariant (must be 0)  (%.0fs)' % (tot, fails, time.time() - T0))
    if os.environ.get('walk', '1') == '1': walkpart(keys)


def walkpart(keys):
    sys.argv = ['x', 'control']
    src = open('workshop/rounds/123/skeptic_mergereplay.py').read()
    g = {'__name__': 'sk'}; exec(compile(src.split('def replay(')[0], 'sk', 'exec'), g)
    from quivermutation import invariants
    r = g['nk'].LinearNakayamaAlgebra.fromClassName('003033030'); pk0 = g['search']._coxeterKeyOrNone(r); x = r
    prev = invs([list(map(int, row)) for row in invariants.integerCartanMatrix(x)])
    print('\n== E-318 walk, invariant changes per step')
    for kind, v in g['PP']('F3 R8 F3 F3 F6'):
        c, lab = g['take'](x, kind, v, pk0)
        cur = invs([list(map(int, row)) for row in invariants.integerCartanMatrix(c)])
        ch = [k for k in keys if prev[k] != cur[k]]
        print('  %s%d  %s %s %s  changed: %s' % (kind, v, 'J!=0' if lab['J'] else 'J=0', 'tp' if lab['tp'] else 'NOT-tp',
              'kept' if lab['kept'] else 'MOVED', ', '.join('%s %s->%s' % (k, prev[k], cur[k]) for k in ch) or 'none'))
        x, prev = c, cur
    print('  end invariants equal to LNA 030330030:', prev == invs(lna(tuple(int(c) for c in '030330030'))), ' (%.0fs)' % (time.time() - T0))


if __name__ == '__main__':
    main()
