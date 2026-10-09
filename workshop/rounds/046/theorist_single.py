"""One placement: pair (P:3)(P+1:3), P=8, plus bystander list given as t:m,...
Usage: theorist_single.py STEPS MARGIN t:m [t:m ...]   (pair at 8,9 always)
Prints reached LNAs, base and min max-overlap, lowered count."""
import sys, time
from quivermutation import lnaMoves as lm, overlap as ov
steps, margin = int(sys.argv[1]), int(sys.argv[2])
by = [tuple(int(x) for x in a.split(':')) for a in sys.argv[3:]]
import os
pair = [tuple(int(x) for x in a.split(':')) for a in os.environ.get('PAIR','8:3,9:3').split(',')]
pat = pair + by
lo = min(s for s, _ in pat); hi = max(s + m - 1 for s, m in pat)
L = hi + margin + 6 + max(0, 8 - lo)
off = max(0, 6 - lo + 0)  # keep >= 6 empty vertices on the left
pat = [(s + off, m) for s, m in pat]; lo += off; hi += off
L = hi + margin + 6
rl = lm.embedPattern(L, pat, 0)
t0 = time.time()
reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
base = ov.maxOverlap(rl)
ovs = [ov.maxOverlap([int(c) for c in n]) for n in reached]
print(f"pattern {pat} L={L} steps={steps} margin={margin} base={base} reached={len(reached)} min={min(ovs) if ovs else '-'} lowered={sum(o<base for o in ovs)} ({time.time()-t0:.0f}s)")
