"""Round 059 (toolsmith), T10 item 1 (E-177 request).  Usage (repo root; needs /tmp/tsm/c<CLS>.pkl, see rounds/058/skeptic_wrongalg.py header):
  CLS=1 timeout 10m .venv/bin/python workshop/rounds/059/toolsmith_dimcheck.py
(1) compare2 (rounds/057/toolsmith_endt2.py, now with dimcheck) on every failing step: dim K Q/I' from crels must equal the child's dims -> assert.
(2) wrong-algebra control: drop each single relation of crels in turn; report dimcheck verdict (expected: rejected iff the dropped relation is not implied by the others)
    and, for comparison, what the old symcheck2 alone said.  Also a spurious-relation control: add the zero-sum 'every path of length 2 between two vertices' is NOT done (not a natural wrong algebra)."""
import pickle, collections, time
exec(open('workshop/rounds/057/toolsmith_endt2.py').read())
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
tab = collections.Counter(); rows = []
for n_, r in enumerate(recs):
    a = r['parentObj']; v = r['v']; c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(a, v))
    TE = TiltEnd(a, v); dims, arrows, layers = endquiver(TE); A, d, carrows, crels = child_info(c)
    cc = collections.Counter((s, h) for (s, h, k) in carrows); par = sorted(m for m in cc.values() if m > 1)
    msg = dimcheck(carrows, crels, d); assert msg is None, (n_, msg)      # (1) the assertion
    tab[('true child', 'dimcheck passes')] += 1
    for i in range(len(crels)):
        new = crels[:i] + crels[i + 1:]; msg = dimcheck(carrows, new, d)
        qd = quotient_dims(carrows, new); big = sum(qd.values()) - sum(d.values())
        t0 = time.time(); old = symcheck2(TE, arrows, 'same', carrows, new); dt = time.time() - t0
        tab[('dropped', 'dim+%d' % big if big else 'dim equal (redundant)', 'dimcheck=' + ('reject' if msg else 'accept'), 'symcheck2 alone=' + ('iso' if old.startswith('iso') else 'NO'))] += 1
        rows.append((n_, v, 'par' + str(par), 'rel', i, len(crels), big, bool(msg), old.split(' (')[0], round(dt, 1)))
    print('step', n_, 'v', v, 'parallel', par, 'nrels', len(crels), 'ok', flush=True)
print('CLS', CLS, 'steps', len(recs))
for k, x in sorted(tab.items(), key=str): print(k, x)
print('rows (step, v, par, _, relidx, nrels, extra dim, dimcheck rejects, symcheck2 alone, s):')
for x in rows: print(x)
