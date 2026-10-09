"""Round 053 (toolsmith), referee items 2 and 3.  Usage (repo root): .venv/bin/python workshop/rounds/053/toolsmith_endt_control.py [gen|wrong]
gen   : per edge of the 13, do the lifted arrows of End(T) generate End(T)?  (span of all products of arrow lifts vs dim Hom(i,j), every i,j)
wrong : per edge, End(T) at (a,v) against OTHER algebras on the same vertex set: (i) the mutation algebra at every other vertex w of a,
        (ii) the 12 other edges' next algebras c', (iii) the 16 failing-step children.  Reports how many pass the dims+arrows filter
        (non-trivial candidates) and how many of those the full relation/Groebner decision rejects; label-preserving, same direction."""
import sys
MODE = sys.argv[1] if len(sys.argv) > 1 else 'gen'
CMODE = MODE
sys.argv = ['x', 'ctl']
_run = open('workshop/rounds/053/toolsmith_endt_run.py').read()
exec(_run.split("if MODE == 'selftest':")[0])
_e = _run[_run.index('def edges13'):_run.index("if MODE == 'perturb':")]
exec(_e)
recs = [r for r in pickle.load(open('/tmp/tsm/c1.pkl', 'rb')) if r['kind'] == 'fail']
E13 = edges13()

if CMODE == 'gen':
    tot = 0
    for name, s_, a, v, c in E13:
        TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); V = TE.V
        lift = {}
        for (i, j, n, t) in arrows: lift.setdefault((i, j), []).append(TE.basis_maps(i, j)[t])
        # S[(i,j)] = echelon of coordinate vectors of products of arrow lifts i ~> j (length >= 1); id for i = j
        S = {(i, j): Echelon() for i in V for j in V}
        cur = {}   # new elements (maps) per (i,j) in the last round
        for (i, j), L in lift.items():
            for f in L:
                if S[(i, j)].add(TE.coords(i, j, f)): cur.setdefault((i, j), []).append(f)
        while cur:
            nxt = {}
            for (i, l), fs in cur.items():
                for (l2, j), gs in lift.items():
                    if l2 != l: continue
                    for f in fs:
                        for g in gs:
                            h = TE.compose(i, l, j, f, g)
                            if h and S[(i, j)].add(TE.coords(i, j, h)): nxt.setdefault((i, j), []).append(h)
            cur = nxt
        ok = all(len(S[(i, j)].rows) == dims[(i, j)] for i in V for j in V if i != j) and all(dims[(i, i)] == 1 for i in V)
        tot += ok
        print(name, s_, 'arrows generate End(T):', ok, '| sum dim End(T) =', sum(dims.values()), '| generated =', len(V) + sum(len(S[k].rows) for k in S if k[0] != k[1]), flush=True)
    print('GEN', tot, 'of', len(E13)); sys.exit()

if CMODE == 'wrong':
    import io, contextlib
    def verdict(a, v, cc):
        if sorted(cc.vertices()) != sorted(a.quiver.nodes): return 'vertexset'
        r = compare(a, v, cc)
        if r['dims'] != 'same': return 'dims'
        if r['arrows'] != 'same': return 'arrows'
        return 'PASS-filter:' + ('iso' if r['rels'].startswith(('all', 'non-monomial relations')) or r.get('sym', '').startswith('iso') else ('undecided' if 'undecided' in r['rels'] + r.get('sym', '') or 'not decided' in r['rels'] else 'REJECT'))
    pool_fail = []
    for r in recs:
        pool_fail.append(r['childObj'])
    T = dict(cand=0, dims=0, arrows=0, filt=0, reject=0, iso=0, undec=0, vset=0)
    for n_, (name, s_, a, v, c) in enumerate(E13):
        cands = []
        for w in sorted(a.quiver.nodes):
            if w == v or not mutation.mutationIsPossibleAtVertex(a, w): continue
            cands.append(('other-vertex w=%d' % w, reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, w))))
        for m_, (n2, s2, a2, v2, c2) in enumerate(E13):
            if m_ != n_: cands.append(('edge %s %s' % (n2, s2), c2))
        for k, f in enumerate(pool_fail): cands.append(('fail-child %d' % k, f))
        cnt = dict(dims=0, arrows=0, filt_iso=0, filt_reject=0, filt_undec=0, vset=0)
        for lab, cc in cands:
            vd = verdict(a, v, cc)
            if vd == 'vertexset': cnt['vset'] += 1
            elif vd in ('dims', 'arrows'): cnt[vd] += 1
            elif vd.endswith('iso'): cnt['filt_iso'] += 1; print('   ISO with', lab)
            elif vd.endswith('REJECT'): cnt['filt_reject'] += 1
            else: cnt['filt_undec'] += 1
        print(name, s_, 'candidates', len(cands), cnt, flush=True)
        for k_, x in cnt.items(): T[{'filt_iso':'iso','filt_reject':'reject','filt_undec':'undec'}.get(k_, k_)] = T.get({'filt_iso':'iso','filt_reject':'reject','filt_undec':'undec'}.get(k_, k_), 0) + x
        T['cand'] += len(cands)
    print('WRONG', T)
