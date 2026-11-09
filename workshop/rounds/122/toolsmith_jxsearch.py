"""Round 122 (toolsmith): T10 / H-015, H-032.  Search guarded balls of NO-capable key groups for a gate-admitted, key-kept J != 0 step whose
child is NOT joined to its parent by J = 0 steps, and test each such step with a Cartan-level invariant; planted crossing edge as power control.

Groups: every LNA key group at length n whose LNAs fall into >= 2 classes of INV = (SNF of S = C + C^T, signature of S, q_S histogram on the
discriminant group when det S != 0).  Different INV => not derived equivalent (E-306, E-308 chain), although the key is equal.
Walk: guarded BFS from ROOT and from its dual, library semantics as `maverick_jtilt.py` (round 121, skeptic): gate
`mutationIsPossibleAtVertex`, illegal-relation drop, `reducePathAlgebra`, key guard (child kept iff its key == root key), cyclic parents not
expanded, dedup by `canonicalKey`.  EVERY key-kept child is expanded, J != 0 ones included (as the library walk does).
Step label: J != 0 iff `perI` (round 033) is non-empty.  (`perI` and `tiltingPlus` of round 001 compute the same kernel, so J = 0 <=> tilting.)
J = 0 component: union-find over canonical keys with every J = 0 key-kept step as an undirected edge, built within the ball only.  For each
key-kept J != 0 step: (a) INV(parent) == INV(child)?  (b) are parent and child in one union-find class at the end of the ball ("joined in
ball") or not ("open": not joined within the ball, the case T10 needs).  A step with INV different is a key-kept step between derived classes.

Controls (CTRL=planted): after the depth-1 layer one PLANTED J != 0 edge is added from the first expanded depth-1 node to an LNA of another INV
class of the group (key equal by construction, so the guard admits it), and one planted in-component edge from that node to one of its own J = 0
children.  The crossing edge must come out (kept, INV differs, open); the in-component one (kept, INV equal, joined).

usage (repo root):
  .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxsearch.py --plan N [N ...]      # groups, classes, depth-2 probe per root, projection
  ROOT=... DEPTH=d [CTRL=planted] timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxsearch.py run N
  ROOTS=all|k DEPTH=d timeout 10m ... run N   # every (or the first k per class) LNA of every NO-capable group at n, one ball each"""
import sys, os, time, copy, collections
from fractions import Fraction
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))          # perI, nk, mutation, procedure, reduction, ap, search, fingerprint, pathAlgebra, nx
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct, invariants

KINDS = os.environ.get('KINDS', 'FR')   # F = right mutation (the library walk), R = left (mutate the dual, dualise back; E-212 `moves`)
_inv = {}
def INV(C):
    C = [list(map(int, r)) for r in C]; k = tuple(map(tuple, C))
    if k in _inv: return _inv[k]
    n = len(C); S = [[C[i][j] + C[j][i] for j in range(n)] for i in range(n)]
    M = sp.Matrix(S); snf = smith_normal_form(M, domain=sp.ZZ); snf = tuple(sorted(abs(int(snf[i, i])) for i in range(n)))
    ev = np.linalg.eigvalsh(np.array(S, dtype=float)); sg = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
    hist = None
    d = int(M.det())
    if d != 0 and abs(d) <= 4096:
        Mi = M.inv(); cols = [[Fraction(int(Mi[i, j].p), int(Mi[i, j].q)) for i in range(n)] for j in range(n)]
        kf = lambda v: tuple(x % 1 for x in v); seen = {kf([Fraction(0)] * n)}; fr = list(seen)
        while fr:
            nf = []
            for g in fr:
                for c in cols:
                    h = kf([a + b for a, b in zip(g, c)])
                    if h not in seen: seen.add(h); nf.append(h)
            fr = nf
        hist = tuple(sorted(collections.Counter(sum(y[i] * S[i][j] * y[j] for i in range(n) for j in range(n)) % 2 for y in seen).items()))
    _inv[k] = (snf, sg, hist); return _inv[k]

def groups(N):
    """NO-capable key groups at length N: list of (key, {INV: [className, ...]})."""
    out = []
    for key, mem in ct.lnaKeyIndex(N).items():
        if len(mem) < 2: continue
        cl = collections.defaultdict(list)
        for rl in mem: cl[INV(ct.lnaCartanMatrix(N, rl))].append(ct.className(rl))
        if len(cl) > 1: out.append((key, dict(cl)))
    return out

class UF:
    def __init__(s): s.p = {}
    def f(s, x):
        s.p.setdefault(x, x)
        while s.p[x] != x: s.p[x] = s.p[s.p[x]]; x = s.p[x]
        return x
    def u(s, a, b): s.p[s.f(a)] = s.f(b)

def lna_name(alg, N):
    q = alg.quiver
    if q.number_of_edges() != N - 1 or len(q.nodes) != N or not nx.is_directed_acyclic_graph(q) or nx.dag_longest_path_length(q) != N - 1: return None
    L = search._renumberedLine(alg); rl = [0] * (N - 2)
    for rel in L.rels:
        if len(rel) != 1: return None
        p = rel[0]; rl[p[0] - 1] = len(p) - 1
    return ct.className(rl)

def ball(N, root, D, label, ctrl='', verbose=True):
    """returns dict of tallies and the list of key-kept J != 0 steps."""
    t0 = time.time(); r = nk.LinearNakayamaAlgebra.fromClassName(root)
    pk0 = search._coxeterKeyOrNone(r); inv0 = INV(invariants.integerCartanMatrix(r))
    seen = set(); frontier = []; tab = collections.Counter(); uf = UF(); jx = []; lnas = {}
    def add(alg, k):
        if k is not None:
            if k in seen: return
            seen.add(k); nm = lna_name(alg, N)
            if nm: lnas[k] = nm
        frontier.append((alg, k))
    for a in (copy.deepcopy(r), pathAlgebra.dualPathAlgebra(r)):
        k = fingerprint.canonicalKey(a); add(a, k); uf.f(k)
    uf.u(fingerprint.canonicalKey(r), fingerprint.canonicalKey(pathAlgebra.dualPathAlgebra(r)))   # root ~ its dual (opposite algebra)
    planted_done = False
    for level in range(D):
        cur, frontier[:] = frontier[:], []
        for alg, pk in cur:
            if list(nx.simple_cycles(alg.quiver)): tab['cyclic parent'] += 1; continue
            pinv = None; j0kids = []
            for kind, v in [(kd, v) for kd in KINDS for v in sorted(alg.quiver.nodes)]:
                a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
                if not mutation.mutationIsPossibleAtVertex(a, v): continue
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): tab['illegal'] += 1; continue
                ch = reduction.reducePathAlgebra(raw)
                if kind == 'R': ch = pathAlgebra.dualPathAlgebra(ch)
                ck = search._coxeterKeyOrNone(ch)
                if ck is None: tab['cyclic-child'] += 1; continue
                J = bool(perI(a, v)); kept = ck == pk0
                tab[('J!=0' if J else 'J=0', 'kept' if kept else 'moved')] += 1
                if not kept: continue
                k = fingerprint.canonicalKey(ch)
                if not J:
                    if pk is not None and k is not None: uf.u(pk, k)
                    j0kids.append((ch, k))
                else:
                    if pinv is None: pinv = INV(invariants.integerCartanMatrix(alg))
                    jx.append(dict(level=level + 1, v=kind + str(v), pk=pk, ck=k, inv_eq=INV(invariants.integerCartanMatrix(ch)) == pinv,
                                   inv_root=INV(invariants.integerCartanMatrix(ch)) == inv0, planted=None, child=ch))
                add(ch, k)
            if ctrl == 'planted' and level == 1 and not planted_done and j0kids:
                planted_done = True
                # crossing: an LNA of another INV class (same key by construction)
                other = label['other']; X = nk.LinearNakayamaAlgebra.fromClassName(other)
                xk = fingerprint.canonicalKey(X); kept = search._coxeterKeyOrNone(X) == pk0
                pinv = INV(invariants.integerCartanMatrix(alg))
                jx.append(dict(level=level + 1, v='planted-cross', pk=pk, ck=xk, inv_eq=INV(invariants.integerCartanMatrix(X)) == pinv,
                               inv_root=INV(invariants.integerCartanMatrix(X)) == inv0, planted='cross', guard_admits=kept))
                if kept: add(X, xk)
                ch, k = j0kids[0]   # in-component: the node's own J = 0 child, relabelled as a J != 0 edge
                jx.append(dict(level=level + 1, v='planted-same', pk=pk, ck=k, inv_eq=INV(invariants.integerCartanMatrix(ch)) == pinv,
                               inv_root=True, planted='same', guard_admits=True, child=ch))
        if verbose: print('  %s level %d distinct %d jx %d %s %.0fs' % (root, level + 1, len(seen), len(jx), dict(tab), time.time() - t0), flush=True)
    EXT = int(os.environ.get('EXT', '3'))
    for s in jx:
        s['joined'] = (s['pk'] is not None and s['ck'] is not None and uf.f(s['pk']) == uf.f(s['ck']))
        s['ext'] = 0 if s['joined'] else None
        if s['joined'] or s['ck'] is None or s['pk'] is None or s['planted'] == 'cross' and False: continue
        # extension: J = 0 key-kept BFS from the child, radius EXT; joined if it meets the parent's in-ball J = 0 class
        if s['planted'] == 'cross': X0 = nk.LinearNakayamaAlgebra.fromClassName(label['other'])
        else: X0 = s['child']
        target = uf.f(s['pk']); fr = [X0]; met = set([s['ck']])
        for e in range(1, EXT + 1):
            nf = []
            for x in fr:
                if list(nx.simple_cycles(x.quiver)): continue
                for kind, v in [(kd, v) for kd in KINDS for v in sorted(x.quiver.nodes)]:
                    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
                    if not mutation.mutationIsPossibleAtVertex(a, v) or perI(a, v): continue
                    raw = mutation.quiverMutationAtVertex(a, v)
                    if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                    ch = reduction.reducePathAlgebra(raw)
                    if kind == 'R': ch = pathAlgebra.dualPathAlgebra(ch)
                    if search._coxeterKeyOrNone(ch) != pk0: continue
                    k = fingerprint.canonicalKey(ch)
                    if k is None or k in met: continue
                    met.add(k)
                    if k in uf.p and uf.f(k) == target: s['ext'] = e; break
                    nf.append(ch)
                if s['ext']: break
            if s['ext']: break
            fr = nf
        s['ext_nodes'] = len(met)
    comp_lnas = collections.defaultdict(set)
    for k, nm in lnas.items(): comp_lnas[uf.f(k)].add(nm)
    for s in jx: s['child_comp_lnas'] = sorted(comp_lnas.get(uf.f(s['ck']), ())) if s['ck'] is not None else None
    return dict(tab=tab, jx=jx, distinct=len(seen), secs=time.time() - t0, lnas=sorted(set(lnas.values())))

def labels(N):
    lab = {}
    for gi, (key, cl) in enumerate(groups(N)):
        ks = list(cl)
        for ci, inv in enumerate(ks):
            for nm in cl[inv]: lab[nm] = dict(group=gi, cls=ci, size=len(cl[inv]), other=cl[ks[(ci + 1) % len(ks)]][0], nclasses=len(ks), gsize=sum(map(len, cl.values())))
    return lab

if __name__ == '__main__':
    if ARGV[1] == '--plan':
        for N in map(int, ARGV[2:]):
            t0 = time.time(); G = groups(N)
            print('n=%d: %d NO-capable key groups (INV-split) [%.0fs]' % (N, len(G), time.time() - t0), flush=True)
            for gi, (key, cl) in enumerate(G):
                print('  group %d: %d LNAs, classes %s' % (gi, sum(map(len, cl.values())), [len(m) for m in cl.values()]), flush=True)
                for inv, mem in cl.items():
                    rt = mem[0]; b = ball(N, rt, 2, None, verbose=False); b1 = ball(N, rt, 1, None, verbose=False)
                    g = b['distinct'] / max(b1['distinct'], 1)
                    print('    class size %d root %s: depth1 %d depth2 %d nodes (%.1fs), growth x%.1f -> depth4 ~%d nodes ~%.0fs, depth5 ~%d ~%.0fs; jx at d<=2: %d' % (
                        len(mem), rt, b1['distinct'], b['distinct'], b['secs'], g, b['distinct'] * g * g, b['secs'] * g * g,
                        b['distinct'] * g ** 3, b['secs'] * g ** 3, len(b['jx'])), flush=True)
        sys.exit()
    N = int(ARGV[2]); D = int(os.environ.get('DEPTH', '3')); ctrl = os.environ.get('CTRL', '')
    lab = labels(N)
    if os.environ.get('ROOT'): roots = os.environ['ROOT'].split(',')
    else:
        sel = os.environ.get('ROOTS', '1'); per = collections.defaultdict(list)
        for nm, l in sorted(lab.items()): per[(l['group'], l['cls'])].append(nm)
        roots = [nm for k in sorted(per) for nm in (per[k] if sel == 'all' else per[k][:int(sel)])]
        if os.environ.get('GROUP'): roots = [r for r in roots if lab[r]['group'] == int(os.environ['GROUP'])]
    T = collections.Counter(); t0 = time.time(); steps = 0
    for rt in roots:
        l = lab[rt]; print('ROOT %s n=%d group %d (%d LNAs) class %d (%d LNAs) depth %d ctrl %s' % (rt, N, l['group'], l['gsize'], l['cls'], l['size'], D, ctrl or '-'), flush=True)
        b = ball(N, rt, D, l, ctrl)
        steps += sum(v for k, v in b['tab'].items() if isinstance(k, tuple))
        for s in b['jx']:
            tag = (s['planted'] or 'real', 'INV eq' if s['inv_eq'] else 'INV DIFFERS', 'joined in ball' if s['joined'] else 'joined at ext %d' % s['ext'] if s['ext'] else 'OPEN')
            T[tag] += 1
            if s['planted'] or not s['joined'] or not s['inv_eq']:
                print('   step level %s v %s planted %s INV(parent)==INV(child) %s INV(child)==INV(root) %s joined-in-ball %s ext-join %s (ext nodes %s) guard admits %s child-component LNAs %s' % (
                    s['level'], s['v'], s['planted'], s['inv_eq'], s['inv_root'], s['joined'], s['ext'], s.get('ext_nodes'), s.get('guard_admits', True), s['child_comp_lnas'][:6]), flush=True)
        print('  ROOT %s: distinct %d, steps %s, key-kept J!=0 %d, LNAs in ball %s, %.0fs' % (rt, b['distinct'], dict(b['tab']), len(b['jx']), b['lnas'][:8], b['secs']), flush=True)
    print('SUMMARY n=%d depth %d roots %d gate-admitted acyclic-child steps %d; key-kept J!=0 steps by (kind, INV, component): %s  %.0fs' % (N, D, len(roots), steps, dict(T), time.time() - t0))
