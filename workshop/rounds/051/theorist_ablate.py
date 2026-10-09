"""T4/H6 ablation: orbit of a two-relation core a,b placed at every offset of length n, reduced walk,
with ALL_MOVES versus floating rules only (anchored rules removed; free=REDUCED, edges/doubles as default).
Verdict = 'inside' iff the walk reaches an almost separate row (stopWhen), else closed/cap.
usage: theorist_ablate.py N A B [limit]"""
import sys
from quivermutation import lnaMoves as lm, freeMoves as fm, overlap
n = int(sys.argv[1]); a = int(sys.argv[2]); b = int(sys.argv[3]); lim = int(sys.argv[4]) if len(sys.argv) > 4 else 40000
for o in range(0, n - 3):
    rl = [0]*(n-2)
    if o + 1 >= len(rl): break
    rl[o] = a; rl[o+1] = b
    if not lm.isAdmissible(n, rl): print(o, 'not admissible'); continue
    out = []
    for name, rules in (('all', None), ('floating', lm.VERIFIED_MOVES), ('norules', [])):
        w = fm.orbitReport(n, rl, rules=rules, free=fm.REDUCED, limit=lim)
        ins = any(overlap.isAlmostSeparate(n, list(r)) for r in w.rows)
        out.append((name, w.stoppedBy, len(w.rows), 'inside' if ins else 'outside', hash(frozenset(w.rows)) % 100000))
    print(o, out, 'SAME' if out[0][1:] == out[1][1:] else 'DIFF', flush=True)
