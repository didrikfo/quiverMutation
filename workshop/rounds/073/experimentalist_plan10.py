"""S-1 / round 073: size an exact classifyLength(N) (default N = 10) without running it.
--plan (default): seed from the quipu theorem + move orbits (the cheap, deterministic first step of
  classifyLength) into a throw-away table, count rows still unplaced, then time `samples` depth-first
  searches from unplaced rows at the depth classifyLength would use (decay: max(depth - floor(log10(i+1)), 2)).
--run: classifyLength with --budget-hours (stops between classes, exit code 2 when the budget is spent,
  0 when finished); resumes from --out. Same as `classify.py N --resume --budget-hours H` plus a plan.
usage: .venv/bin/python workshop/rounds/073/experimentalist_plan10.py [--n 10] [--samples 6] [--sample-seconds 120]
       .venv/bin/python workshop/rounds/073/experimentalist_plan10.py --run --budget-hours 9 --out logs/A_10_exact.csv"""
import argparse, math, sys, time, tempfile, os, collections
from quivermutation import classification as cl, mutationClassTable as mct, lines, nakayama, search

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=10); ap.add_argument('--depth', type=int, default=6)
    ap.add_argument('--samples', type=int, default=6); ap.add_argument('--sample-seconds', type=float, default=120)
    ap.add_argument('--run', action='store_true'); ap.add_argument('--budget-hours', type=float, default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.run:
        out = a.out or 'A_%d_mutation_classes.csv' % a.n
        t, rep = cl.classifyLength(a.n, a.depth, 6, out, printOutput=False, resume=True,
                                   budgetSeconds=None if a.budget_hours is None else a.budget_hours * 3600)
        print('rows', len(t), 'unplaced', rep['unplacedRows'], 'classes', len(t.classNames()),
              'candidates', len(rep['candidate']), 'stoppedEarly', rep['stoppedEarly'])
        sys.exit(2 if rep['stoppedEarly'] else 0)
    n = a.n
    t0 = time.time()
    tmp = os.path.join(tempfile.mkdtemp(), 'plan.csv')
    table = mct.MutationClassTable.forLength(n, [lines.relSetToString(r) for r in lines.generateAllPossibleLineRelations(n)])
    seeded = cl.seedTableFromQuipuTheorem(table, n, printOutput=False)
    un = [i for i, r in enumerate(table.rows()) if not r[1]]
    print('n=%d rows %d; placed by theorem+moves %d (%.1f%%) in %d classes; unplaced %d; seed step %.0f s'
          % (n, len(table), len(table) - len(un), 100 * (len(table) - len(un)) / len(table), len(table.classNames()), len(un), time.time() - t0))
    dec = collections.Counter(int(max(a.depth - math.floor(math.log10(i + 1)), 2)) for i in un)
    print('unplaced rows by search depth (upper bound, decay by row index):', dict(sorted(dec.items())))
    # time sample searches, evenly spaced over the unplaced rows, each with a wall-clock cap via depth-first on a thread-less check
    import signal
    class TO(Exception): pass
    def h(s, f): raise TO()
    signal.signal(signal.SIGALRM, h)
    step = max(1, len(un) // a.samples)
    for i in un[::step][:a.samples]:
        d = int(max(a.depth - math.floor(math.log10(i + 1)), 2))
        al = nakayama.LinearNakayamaAlgebra.fromRelationString(n, table.rows()[i][0])
        out = []; t1 = time.time()
        signal.alarm(int(a.sample_seconds))
        try:
            search.mutationSearchDepthFirst(al, d, [], 'x', printOutput=False, collected=out)
            r = 'done'
        except TO:
            r = 'CAPPED'
        signal.alarm(0)
        print('row %5d depth %d: %-6s %6.1f s, %d LNAs reached' % (i + 1, d, r, time.time() - t1, len(out)))
main()
