"""Scholar r051: try to find a 7-step path 05040330 -> 33460000 (n=10, E-033 link A), split by first mutation, then replay with tiltingPlus.
Usage: timeout 9m .venv/bin/python workshop/rounds/051/scholar_apath.py <startrow-as-digits> <depth>"""
import sys, copy, time, multiprocessing as mp
sys.path.insert(0, '.')
from quivermutation import lnaMoves as lm, nakayama as nk, mutation, procedure, reduction, search, fingerprint
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'))
N = 10
def work(a):
    row, depth, si, v = a
    rs = nk.LinearNakayamaAlgebra(N, row).relationString()
    start = list(search.memberAndItsDual(N, rs))[si]
    if not mutation.mutationIsPossibleAtVertex(start, v): return []
    alg0 = copy.deepcopy(start)
    child = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(start), v))
    col = []; vis = fingerprint.Visited()
    search.mutationSearchDepthFirst(child, depth - 1, [], 'x', printOutput=False, collected=col, visited=vis)
    out = []
    for alg, path, _ in col:
        r = lm.asRelLengths(lm._copy(alg), N)
        if r is not None: out.append((si, [v] + list(path), tuple(r)))
    return out
if __name__ == '__main__':
    row = [int(c) for c in sys.argv[1]]; depth = int(sys.argv[2])
    tasks = [(row, depth, si, v) for si in (0, 1) for v in range(1, N + 1)]
    t = time.time()
    with mp.Pool(4) as p:
        for res in p.imap_unordered(work, tasks):
            for r in res: print(r, flush=True)
            print('# done', round(time.time() - t), flush=True)
