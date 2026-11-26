"""Round 139 (toolsmith), T10 / F-057 (a): Hom check of the nine J = 0 steps of the E-332 10-step walk.
Pre-registration: workshop/rounds/139/toolsmith_prereg.txt (written before the first run).

Per step (parent A, vertex k; R steps as F on the dual) three verdicts:
  LIB = mutation.isTiltingAtVertex(A, k)
  HOM = Hom complex of T = Y + sum_{i != k} P_i, Y = [P_k -> sum_beta P_t(beta)], by the round-129 theorist engine
        (workshop/rounds/129/theorist_jsilt.py, imported unchanged): H[m], m = -2, -1, 1, 2, all 0; End^0(T) Cartan vs child
  OWN = new here: own path enumeration, ideal = span of all p r q (exact rationals, no library ideal basis), then
        Hom(T,T[-1]) = sum_{i != k} dim ker(e_i A e_k -> sum_beta e_i A e_t(beta)) + Hom(Y,Y[-1]) (= sum_beta dim e_t A e_k),
        Hom(T,T[1]) = 0 iff the arrows out of k left-approximate (every e_k A e_j, j != k, in sum_beta beta . e_t A e_j).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/139/toolsmith_homj0.py [walks|controls|e316|all]
"""
import sys, importlib.util, collections, time, warnings
from fractions import Fraction
sys.path.insert(0, '.')
import networkx as nx
from quivermutation import arrowPaths as ap, procedure, pathAlgebra, mutation, reduction, search, nakayama as nk
warnings.simplefilter('ignore')

_spec = importlib.util.spec_from_file_location('jsilt', 'workshop/rounds/129/theorist_jsilt.py')
TH = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(TH)

# ------------------------------------------------------------------ own code
def rank(rows):
    piv = {}; r = 0
    for row in rows:
        row = {k: Fraction(v) for k, v in row.items() if v}
        while row:
            h = min(row, key=repr)
            if h not in piv: piv[h] = {k: v / row[h] for k, v in row.items()}; r += 1; break
            f, pv = row[h], piv[h]
            row = {k: row.get(k, 0) - f * pv.get(k, 0) for k in set(row) | set(pv)}
            row = {k: v for k, v in row.items() if v}
    return r

class Own:
    def __init__(s, A):
        s.Q = A.quiver
        assert nx.is_directed_acyclic_graph(s.Q), 'acyclic only'
        s.E = sorted(s.Q.edges(keys=True)); s.V = sorted(s.Q.nodes)
        s.P = collections.defaultdict(list)            # (i, j) -> paths (tuples of edges); () only for i == j
        for v in s.V: s.P[(v, v)].append(())
        stack = [(e,) for e in s.E]
        while stack:
            p = stack.pop(); s.P[(p[0][0], p[-1][1])].append(p)
            stack += [p + (e,) for e in s.E if e[0] == p[-1][1]]
        s.I = collections.defaultdict(list)            # (i, j) -> spanning vectors of e_i I e_j
        for r in procedure.relationsFrom(A):
            terms = [(tuple(tuple(a) for a in p), Fraction(c)) for p, c in r.items() if c]
            a = {p[0][0] for p, _ in terms}; b = {p[-1][1] for p, _ in terms}
            assert len(a) == 1 and len(b) == 1, 'inhomogeneous relation'
            a, b = a.pop(), b.pop()
            for i in s.V:
                for p in s.P[(i, a)]:
                    for j in s.V:
                        for q in s.P[(b, j)]:
                            vec = collections.Counter()
                            for w, c in terms: vec[p + w + q] += c
                            s.I[(i, j)].append(dict(vec))
    def dim(s, i, j): return len(s.P[(i, j)]) - rank(s.I[(i, j)])

def own_verdict(A, k):
    O = Own(A); outs = [e for e in O.E if e[0] == k]
    Hm1 = 0
    for i in O.V:                                    # Hom(P_i, Y[-1]) = ker of right mult by the arrows out of k
        if i == k or not O.P[(i, k)]: continue
        Ib = [{(e, p): c for p, c in v.items()} for e in outs for v in O.I[(i, e[1])]]
        img = [{(e, x + (e,)): 1 for e in outs} for x in O.P[(i, k)]]
        rk_mod = rank(img + Ib) - rank(Ib)           # rank of the induced map from paths into the quotients
        Hm1 += len(O.P[(i, k)]) - rk_mod - rank(O.I[(i, k)])
    HYY = sum(O.dim(e[1], k) for e in outs)         # Hom(Y, Y[-1]) = Hom(sum P_t, P_k); 0 on an acyclic quiver
    appr = True
    for j in O.V:
        if j == k or not O.P[(k, j)]: continue
        sub = [{(e,) + q: 1} for e in outs for q in O.P[(e[1], j)]]
        if rank(sub + O.I[(k, j)]) < rank([{p: 1} for p in O.P[(k, j)]] + O.I[(k, j)]): appr = False
    return dict(own_Hm1=Hm1 + HYY, own_HYY=HYY, own_appr=appr, own_tilt=(Hm1 + HYY == 0 and appr))

# ------------------------------------------------------------------ steps
def child_of(A, k):
    assert mutation.mutationIsPossibleAtVertex(A, k), 'gate refuses %d' % k
    raw = mutation.quiverMutationAtVertex(A, k)
    assert not any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)), 'illegal child relation'
    return reduction.reducePathAlgebra(raw)

TALLY = collections.Counter(); DIS = []

def check(label, A, k, child, cls):
    t = TH.test(A, k, child)
    o = own_verdict(A, k)
    illegal = any(ap.isIllegalRelation(A.quiver, r) for r in procedure.relationsFrom(A))
    hom = t['tilting']; lib = t['lib_tilt']; own = o['own_tilt']
    agree = lib == hom == own and o['own_Hm1'] == t['H'][-1] and o['own_appr'] == t['appr']
    if not agree: DIS.append(label)
    TALLY[(cls, 'J=%d' % t['J'], 'LIB %s' % lib, 'HOM %s' % hom, 'OWN %s' % own, 'agree %s' % agree,
           'End0~child %s' % t['End0_eq_child'], 'silting %s' % t['silting'])] += 1
    print('%-44s J %d | LIB %-5s | HOM H[-2,-1,1,2] %s %s %s %s tilt %-5s | OWN H[-1] %d (YY %d) appr %-5s tilt %-5s | End0~child %-5s illegalRel %s%s' % (
        label, t['J'], lib, t['H'][-2], t['H'][-1], t['H'][1], t['H'][2], hom, o['own_Hm1'], o['own_HYY'], o['own_appr'], own,
        t['End0_eq_child'], illegal, '' if agree else '   <-- DISAGREEMENT'), flush=True)
    return t, o

def replay(root, side, moves, cls, label, check_all=True, last_cls=None):
    X = nk.LinearNakayamaAlgebra.fromClassName(root)
    if side == 'dual': X = pathAlgebra.dualPathAlgebra(X)
    out = []
    for d, mv in enumerate(moves, 1):
        k = int(mv[1:]); A = X if mv[0] == 'F' else pathAlgebra.dualPathAlgebra(X)
        ch = child_of(A, k)
        if check_all or d == len(moves):
            c = last_cls if (d == len(moves) and last_cls) else cls
            out.append(check('%s %s(%s) d%d %s' % (label, side, root, d, mv), A, k, ch, c))
        X = ch if mv[0] == 'F' else pathAlgebra.dualPathAlgebra(ch)
    return X, out

def walks():
    print('== target: E-332 / F-057 (a) 10-step walk, and the second root (5 J = 0 + 1 J != 0)')
    b = 'F6 F11 F7 F11 F7 F7 F3 F3 F4 F7'.split()
    X = nk.LinearNakayamaAlgebra.fromClassName('303305000'); X = pathAlgebra.dualPathAlgebra(X)
    for d, mv in enumerate(b, 1):
        k = int(mv[1:]); ch = child_of(X, k)
        check('walk b dual(303305000) d%-2d %s' % (d, mv), X, k, ch, 'C1 control J!=0 (walk b step 6)' if d == 6 else 'TARGET walk b J=0')
        X = ch
    a = 'F1 F1 F4 F3 F1 F8'.split()
    X = nk.LinearNakayamaAlgebra.fromClassName('303500030')
    for d, mv in enumerate(a, 1):
        k = int(mv[1:]); ch = child_of(X, k)
        check('walk a 303500030 d%-2d %s' % (d, mv), X, k, ch, 'C1 control J!=0 (walk a step 6)' if d == 6 else 'secondary walk a J=0')
        X = ch

def controls():
    print('== C3 / YES control (hand-built, theorist_jsilt.build)')
    sq = [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
    E080 = TH.build(sq, [{(0, 2, 4): 1, (1, 3, 4): -1}]); SQ = TH.build(sq, [{(0, 2): 1, (1, 3): -1}])
    check('C3 E-080 abde=acde @4 F', E080, 4, TH.child_of(E080, 4), 'C3 control E-080')
    check('YES abd=acd @4 F', SQ, 4, TH.child_of(SQ, 4), 'YES control abd=acd')

E316 = [('E-233', '3060000', 'F1 R9 R8 R8', 'R8'),
        ('E-234 1275', '3060000', 'F1 F4 R9 R8 R8', 'R8'), ('E-234 1433', '3060000', 'F1 F5 R9 R8 R8', 'R8'),
        ('E-234 1746', '3060000', 'F1 R5 R9 R8 R8', 'R8'), ('E-234 1818', '3060000', 'F1 R6 R9 R8 R8', 'R8'),
        ('E-234 1915', '3060000', 'F1 R8 R9 F8 R9', 'R9'), ('E-234 1932', '3060000', 'F1 R8 R9 R9 R8', 'R8'),
        ('E-234 1970', '3060000', 'F1 R9 F1 R8 R8', 'R8'), ('E-234 2003', '3060000', 'F1 R9 F3 R8 R8', 'R8'),
        ('E-234 4551', '3060000', 'R8 R9 R9 F1 R8', 'R8'), ('E-234 4722', '3060000', 'R9 R8 R8 F1 F8', 'R1'),
        ('E-239 2026', '3304000', 'F2 R8 R9 R7 R7', 'R7'),
        ('E-244 4606', '3304000', 'F2 F1 R8 R9 R7 R7', 'R7'), ('E-244 8063', '3304000', 'F2 R8 R7 R9 F7 R9', 'R9'),
        ('E-244 8170', '3304000', 'F2 R8 R9 F2 R7 R7', 'R7'), ('E-244 8196', '3304000', 'F2 R8 R9 F4 R7 R7', 'R7'),
        ('E-244 8427', '3304000', 'F2 R9 R8 R7 R9 R7', 'R7'), ('E-244 11123', '3304000', 'R8 R9 R7 R7 F2 F7', 'R2'),
        ('E-254', '2400230', 'F1 F1 F2 F3 R9 R7 R9 F3', 'F3')]   # copied from workshop/rounds/121/skeptic_jcount.py replay

def e316():
    print('== C2: E-316 19 n = 9 key-kept J != 0 steps (path moves checked too: they are the J = 0 side of the joins)')
    for rec, root, path, step in E316:
        try: replay(root, 'root', path.split() + [step], 'E-316 path move (J=0 expected)', rec, True, 'C2 control E-316 step')
        except AssertionError as e: print('%-12s replay broken: %s' % (rec, e)); TALLY[('C2 replay broken',)] += 1

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'; t0 = time.time()
    if mode in ('walks', 'all'): walks()
    if mode in ('controls', 'all'): controls()
    if mode in ('e316', 'all'): e316()
    print('TALLY (class, J, LIB, HOM, OWN, all agree, End0 ~ child Cartan, silting) -> count')
    for k, v in sorted(TALLY.items(), key=str): print('  ', k, v)
    print('DISAGREEMENTS:', len(DIS), DIS)
    print('time %.0fs' % (time.time() - t0))
