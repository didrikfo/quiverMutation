"""Round 111 toolsmith: NO control for toolsmith_c1_13_exact.congruent (no pari): E-269's n = 12 group of 140, cross pair (certified inequivalent, E-259) must say NO, within-quipu YES.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/111/toolsmith_c1_13_ctl12.py"""
import sys; sys.path.insert(0, 'workshop/rounds/111'); sys.path.insert(0, '.')
import toolsmith_c1_13_exact as E, numpy as np, itertools, collections
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
E.N = 12
C = lambda rl: np.array(ct.lnaCartanMatrix(12, rl)).astype(np.int64)
a0 = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); keys = ct.lnaKeys(12); kk = keys[a0]
grp = [rl for rl, k in keys.items() if k == kk]; q = collections.defaultdict(list)
for rl in grp:
    r = quipuForAlmostSeparateLNA(12, list(rl))
    if r is not None: q[r].append(tuple(rl))
print('group', len(grp), 'quipus', {k: len(v) for k, v in q.items()})
reps = [(r, rl) for r in q for rl in q[r][:2]]
for (r1, a), (r2, b) in itertools.combinations(reps, 2):
    d = E.congruent(C(a), C(b), 0.5); print('within' if r1 == r2 else 'ACROSS', a, b, d, flush=True)
