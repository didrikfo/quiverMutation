"""Round 129 (skeptic), T10 blocking check for E-323, part (1): Hom-level test of the F3 step at P = 0040033030.F3 R9 F3.
Hom engine = rounds/091/skeptic_convention.py with the vector-valued (non-monomial) patch of rounds/118/skeptic_conv3.py, mod 32003 --
the engine of E-327, independent of the library's gate, of perI (J) and of isTiltingAtVertex.  For a step (algebra a, vertex v; an R step
is F on the dual, as in toolsmith_jxhit.py) T = (+)_{u != v} P_u (+) [P_v -> (+)_{v->j} P_j] (the 'R' complex of E-327 = F side), and
  tilting(Hom)  <=>  Hom(T, T[m]) = 0 for m = -2, -1, 1, 2   (generation is automatic for this complex; two-term, so |m| >= 2 is a check).
Also: dim Hom(T_x, T_y) (m = 0) against the library child's Cartan matrix (equal up to transpose expected exactly when tilting).
Cases: the E-323 step and its dual mirror (pre-registered NO); the three J = 0 steps of the path F3 R9 F3 from the root (pre-registered YES);
every gate-admitted F/R step at P and at the root (J vs Hom agreement); E-080's algebra abde = acde at d (NO control) and abd = acd (YES).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/129/skeptic_hom323.py"""
import sys, copy, time
MYARGV = sys.argv[1:]
sys.path.insert(0, '.'); sys.argv = ['x', 'none']
_s = open('workshop/rounds/122/toolsmith_jxsearch.py').read(); exec(compile(_s.split("if __name__ == '__main__':")[0], 'jx', 'exec'))
_h = open('workshop/rounds/122/toolsmith_jxhit.py').read(); exec(compile(_h[_h.index('def step('):_h.index('lab = labels(N)')], 'jxhit', 'exec'))
sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
import skeptic_conv3 as c3            # patches skeptic_convention's Dmat / comp_classes to the non-monomial engine
import skeptic_convention as sc
P = sc.P

def toAlgG(alg):
    V = sorted(alg.quiver.nodes); ix = {v: t for t, v in enumerate(V)}
    E = list(alg.quiver.edges(keys=True)); ex = {e: t for t, e in enumerate(E)}
    rels = [{tuple(ex[e] for e in p): int(c) % P for p, c in r.items()} for r in procedure.relationsFrom(alg)]
    assert all(c == int(c) for r in procedure.relationsFrom(alg) for c in r.values())
    return c3.AlgG(len(V), [(ix[a], ix[b]) for a, b, _ in E], rels, True), ix

def homT(A, i):
    summ = [sc.stalk(u) if u != i else sc.mutate_cx(A, i, 'R') for u in range(A.n)]
    deg = {}
    for m in (-2, -1, 1, 2):
        d = sum(sc.cx_hom(A, summ[x], summ[y], m).dim for x in range(A.n) for y in range(A.n)
                if x == i or y == i)            # stalk-to-stalk Hom in degree != 0 vanishes
        if d: deg[m] = d
    M = [[sc.cx_hom(A, summ[x], summ[y], 0).dim for y in range(A.n)] for x in range(A.n)]
    return deg, M

def check(alg, mv, label):
    kind, v = mv[0], int(mv[1:]); a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
    gate = mutation.mutationIsPossibleAtVertex(a, v)
    if not gate: return None
    child, J = step(alg, mv)
    A, ix = toAlgG(a); deg, M = homT(A, ix[v])
    ca = a if kind == 'F' else None
    chF = child if kind == 'F' else pathAlgebra.dualPathAlgebra(child)     # child on the side of a
    C = invariants.integerCartanMatrix(chF) if nx.is_directed_acyclic_graph(chF.quiver) else None
    Ct = None if C is None else [list(r) for r in zip(*C)]
    cart = None if C is None else (M == C or M == Ct)
    res = dict(label=label, mv=mv, J=J, inj=mutation.isTiltingAtVertex(a, v), hom=not deg, deg=deg, cartan=cart,
               key=search._coxeterKeyOrNone(child) == search._coxeterKeyOrNone(alg))
    print('  %-28s %-4s J!=0 %-5s isTiltingAtVertex %-5s Hom-tilting %-5s Hom(T,T[m])!=0 at %-10s End(T)=child Cartan %-5s key kept %s'
          % (label, mv, J, res['inj'], res['hom'], deg, cart, res['key']), flush=True)
    return res

def e080():
    a, b, c, d, e = range(5); arrows = [(a, b), (a, c), (b, d), (c, d), (d, e)]
    for name, rel, i in (('E-080 abde=acde at d (NO)', {(0, 2, 4): 1, (1, 3, 4): P - 1}, d), ('abd=acd at d (YES)', {(0, 2): 1, (1, 3): P - 1}, d)):
        A = c3.AlgG(5, arrows, [rel], True); deg, M = homT(A, i)
        print('  %-28s Hom-tilting %s, Hom(T,T[m])!=0 at %s' % (name, not deg, deg))

if __name__ == '__main__':
    t0 = time.time(); N = 12
    print('controls (E-327 shapes):'); e080()
    r = nk.LinearNakayamaAlgebra.fromClassName('0040033030')
    print('E-323 path from 0040033030 (J = 0 steps, pre-registered Hom YES), then the step F3 (pre-registered Hom NO):')
    x = copy.deepcopy(r); out = []
    for mv in ('F3', 'R9', 'F3'):
        out.append(check(x, mv, 'path step')); x, _ = step(x, mv)
    Ppar = x; out.append(check(Ppar, 'F3', 'E-323 step at P'))
    d = pathAlgebra.dualPathAlgebra(r)
    for mv in ('F9', 'R3', 'R3'): d, _ = step(d, mv)
    out.append(check(d, 'R3', 'E-323 mirror (dual side)'))
    print('relations of P:', procedure.relationsFrom(Ppar))
    for nm, alg in (('root 0040033030', r), ('P', Ppar)):
        print('every gate-admitted step at', nm)
        for kind in 'FR':
            for v in range(1, N + 1):
                res = check(alg, kind + str(v), nm)
                if res: out.append(res)
    bad = [o for o in out if o and (o['J'] == o['hom'] or o['inj'] != o['hom'] or (o['hom'] and o['cartan'] is False))]
    print('cases %d; disagreements (J!=0 vs Hom, isTiltingAtVertex vs Hom, Hom-tilting but End(T) != child Cartan): %d %s' % (len(out), len(bad), [(o['label'], o['mv']) for o in bad]))
    print('%.0fs' % (time.time() - t0))
