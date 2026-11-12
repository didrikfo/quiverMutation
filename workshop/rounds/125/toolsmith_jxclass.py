"""Round 125 (toolsmith): T10 classifier.  For every gate-admitted key-kept J != 0 step in the guarded F/R ball of each root (round-122
`toolsmith_jxsearch.ball`, same semantics: gate, illegal-relation drop, reducePathAlgebra, key guard, every key-kept child expanded), record
  J      : perI non-empty (as the walk labels it),
  tilt   : `tiltingPlus` (round 001, Ladkani Prop 2.3(c)) on the parent at the vertex (on the dual for an R step); J != 0 <=> not tilt is checked,
  SNF    : SNF(S) tail of parent and child (S = C + C^T), sig: signature of S, qS: q_S histogram (None if det S = 0 or |det S| > 4096),
  cross  : INV(parent) != INV(child)  (=> not a derived equivalence, finite gl.dim both since acyclic: E-318 necessity step).
Crossings are deduplicated up to isomorphism of the (parent, child) pair, also against the dual pair, so a step and its mirror count once.
EXT is forced to 0 (no J = 0 extension), so the "joined" label is not computed here; crossing is decided by INV alone.
CTRL=planted adds round 122's planted crossing edge and planted in-component edge after the depth-1 layer (power control: cross YES / NO).

usage (repo root):
  ROOTS=r1,r2,... DEPTH=d [CTRL=planted] [OUT=file.jsonl] timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxclass.py N
  GROUP=g [SLICE=a:b] DEPTH=d ... toolsmith_jxclass.py N        # roots of key group g (cache order), names sorted, slice a:b
  .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxclass.py --tally file.jsonl [...]   # merge outputs, dedup crossings by iso"""
import sys, os, json, time, collections
sys.path.insert(0, '.'); _A = sys.argv[:]; sys.argv = ['x', 'none']
os.environ['EXT'] = '0'
_s = open('workshop/rounds/122/toolsmith_jxsearch.py').read(); exec(compile(_s.split("if __name__ == '__main__':")[0], 'jx', 'exec'))
_h = open('workshop/rounds/001/scholar_h015.py').read(); exec(compile(_h[_h.index('def tiltingPlus'):_h.index('def cartan(')], 'h015', 'exec'))
sys.argv = _A

def tail(inv): return list(inv[0][-3:])

def classify(N, rt, D, l, ctrl):
    b = ball(N, rt, D, l, ctrl, verbose=False); rows = []
    for s in b['jx']:
        r = dict(root=rt, level=s['level'], v=s['v'], planted=s['planted'], path=' '.join(s.get('path', ())[1:]), side=(s.get('path') or ('-',))[0])
        if s['planted'] == 'cross':
            X = nk.LinearNakayamaAlgebra.fromClassName(l['other']); r.update(child_inv=INV(invariants.integerCartanMatrix(X)))
        if s.get('parent') is not None:
            P, C = s['parent'], s['child']; kind, v = s['v'][0], int(s['v'][1:])
            a = P if kind == 'F' else pathAlgebra.dualPathAlgebra(P)
            r['tilt'] = bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v)); r['J'] = bool(perI(a, v))
            pi, ci = INV(invariants.integerCartanMatrix(P)), INV(invariants.integerCartanMatrix(C))
            r.update(parent_snf=tail(pi), child_snf=tail(ci), parent_sig=pi[1], child_sig=ci[1], parent_qs=pi[2], child_qs=ci[2],
                     snf_changes=pi[0] != ci[0], qs_changes=pi[2] != ci[2], cross=pi != ci,
                     acyclic=bool(nx.is_directed_acyclic_graph(P.quiver) and nx.is_directed_acyclic_graph(C.quiver)),
                     parent_struct=_ser(P), child_struct=_ser(C))
        elif s['planted'] == 'same':
            r.update(cross=not s['inv_eq'])
        else:
            r.update(cross=not s['inv_eq'])
        rows.append(r)
    steps = sum(v for k, v in b['tab'].items() if isinstance(k, tuple))
    head = dict(root=rt, n=N, depth=D, group=l['group'], cls=l['cls'], distinct=b['distinct'], steps=steps,
                tab={'/'.join(k) if isinstance(k, tuple) else k: v for k, v in b['tab'].items()}, secs=round(b['secs'], 1), njx=len(rows))
    return head, rows

def _ser(A):
    return dict(arrows=[list(e) for e in A.quiver.edges(keys=True)], rels=[[[list(map(list, p)), c] for p, c in r.items()] for r in procedure.relationsFrom(A)])

def _deser(d):
    from quivermutation import pathAlgebra as pa
    A = pa.PathAlgebra(); Q = nx.MultiDiGraph()
    for x, y, k in d['arrows']: Q.add_edge(x, y, key=k)
    A.quiver = Q
    ar = [ap.combination({tuple(tuple(a) for a in p): c for p, c in r}) for r in d['rels']]
    A.arrowRels = ar; A.rels = [[list(ap.projectPath(p)) for p in r] for r in ar]
    assert [dict(x) for x in procedure.relationsFrom(A)] == [dict(x) for x in ar]
    return A

def tally(files):
    heads, rows = [], []
    for f in files:
        for line in open(f):
            o = json.loads(line); (heads if 'distinct' in o else rows).append(o)
    print('roots %d (%s), nodes %d, gate-admitted steps %d, key-kept J!=0 steps %d' % (len(heads), sorted({(h['n'], h['depth'], h['group']) for h in heads}),
          sum(h['distinct'] for h in heads), sum(h['steps'] for h in heads), sum(1 for r in rows if not r['planted'])))
    real = [r for r in rows if not r['planted']]
    T = collections.Counter((r['tilt'], r['J'], r['snf_changes'], r['qs_changes'], r['cross']) for r in real)
    print('(tilt, J, SNF changes, qS changes, cross) ->', dict(T))
    print('planted:', collections.Counter((r['planted'], r['cross']) for r in rows if r['planted']))
    reps = []   # iso classes of crossing (parent, child) pairs, up to duality
    for r in real:
        if not r['cross']: continue
        P, C = _deser(r['parent_struct']), _deser(r['child_struct'])
        hit = None
        for i, (P2, C2, rs) in enumerate(reps):
            if (iso(P, P2) and iso(C, C2)) or (iso(pathAlgebra.dualPathAlgebra(P), P2) and iso(pathAlgebra.dualPathAlgebra(C), C2)): hit = i; break
        if hit is None: reps.append((P, C, [r]))
        else: reps[hit][2].append(r)
    print('distinct crossing steps up to iso and duality: %d' % len(reps))
    for P, C, rs in reps:
        print('  SNF %s -> %s, sig %s -> %s; seen %d times from roots %s; e.g. %s %s|%s at level %d' % (rs[0]['parent_snf'], rs[0]['child_snf'],
              rs[0]['parent_sig'], rs[0]['child_sig'], len(rs), sorted({x['root'] for x in rs}), rs[0]['root'], rs[0]['side'], rs[0]['path'] + ' ' + rs[0]['v'], rs[0]['level']))
    roots_with = collections.Counter(r['root'] for r in real if r['cross'])
    print('roots whose ball holds a crossing: %d of %d: %s' % (len(roots_with), len(heads), dict(roots_with)))

if __name__ == '__main__':
    if _A[1] == '--tally': tally(_A[2:]); sys.exit()
    N = int(_A[1]); D = int(os.environ.get('DEPTH', '3')); ctrl = os.environ.get('CTRL', ''); lab = labels(N)
    if os.environ.get('ROOTS'): roots = os.environ['ROOTS'].split(',')
    else:
        g = int(os.environ['GROUP']); roots = sorted(nm for nm, l in lab.items() if l['group'] == g)
        a, b = (os.environ.get('SLICE', ':').split(':') + [''])[:2]; roots = roots[int(a) if a else None:int(b) if b else None]
    out = open(os.environ['OUT'], 'a') if os.environ.get('OUT') else None
    for rt in roots:
        head, rows = classify(N, rt, D, lab[rt], ctrl)
        print('ROOT %s g%d c%d depth %d: distinct %d steps %d key-kept J!=0 %d cross %d (%.0fs)' % (rt, head['group'], head['cls'], D, head['distinct'], head['steps'],
              len(rows), sum(1 for r in rows if r['cross'] and not r['planted']), head['secs']), flush=True)
        for r in rows:
            print('   ', {k: v for k, v in r.items() if not k.endswith('struct')}, flush=True)
        if out:
            out.write(json.dumps(head) + '\n')
            for r in rows: out.write(json.dumps(r, default=list) + '\n')
            out.flush()
