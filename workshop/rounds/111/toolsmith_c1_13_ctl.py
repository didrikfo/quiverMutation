"""Round 111 toolsmith: NO-power probe for toolsmith_c1_13_exact.congruent at n = 7: pairs of key-class-1 Cartan matrices with the same Phi (all share the char poly):
all 89 parents x all 89 children cross pairs (distinct matrices), and C vs C^T.  Reports the tally; any NO would need an independent certificate (none: q_S equal)."""
import sys, pickle, collections
sys.path.insert(0, 'workshop/rounds/111'); ARGV = sys.argv
import toolsmith_c1_13_exact as E, numpy as np
seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
rows = [r for r in sorted(byCart.values(), key=lambda r: r['nexp']) if r['kept'] and not r['tp']]
mats = {}
for r in rows:
    for C in (r['CA'], r['CB']): mats[C.tobytes()] = C
Ms = list(mats.values()); print('distinct Cartan matrices', len(Ms))
t = collections.Counter()
for i, A in enumerate(Ms):
    for B in Ms[i + 1:]:
        try: t[E.congruent(A, B)['congruent']] += 1
        except AssertionError: t['assert'] += 1
print('cross pairs', dict(t))
print('C vs C^T', collections.Counter(E.congruent(A, A.T)['congruent'] for A in Ms))
