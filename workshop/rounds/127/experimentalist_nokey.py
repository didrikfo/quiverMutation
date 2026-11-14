"""Round 127 (response to skeptic, item 4): re-run experimentalist_gap2fix cells with
localMutationSequences' dedup key made faithful on parallel-arrow states (each such state gets
a unique key, so none is pruned), or only count them. Re-created from the referee's probe.
Usage: experimentalist_nokey.py patch|count STEPS MARGIN CAP "s,l:s,l;..." [--end]"""
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import experimentalist_gap2fix as g
from quivermutation import lnaMoves as lm, arrowPaths as ap
orig = lm._stateKey
ctr = itertools.count()
stats = {"par": 0}
def key(p):
    if ap.hasParallelArrows(p.quiver):
        stats["par"] += 1
        return ("par", next(ctr))
    return orig(p)
if sys.argv[1] == "patch": lm._stateKey = key
else:
    def countonly(p):
        if ap.hasParallelArrows(p.quiver): stats["par"] += 1
        return orig(p)
    lm._stateKey = countonly
steps, margin, cap = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
end = "--end" in sys.argv
for c in sys.argv[5].split(";"):
    cell = tuple(tuple(map(int, r.split(","))) for r in c.split(":"))
    stats["par"] = 0
    print(sys.argv[1], *g.run((cell, steps, margin, cap, end)), "parallel_states", stats["par"], flush=True)
