"""Parallel (first-move-split) version of lnaMoves.localMutationSequences for T7 k = 5 cells.
Each worker takes one first signed mutation and runs the library's search below it with its own `seen` table;
the union of reached LNA names equals the serial result (seen-sharing across first moves is only a speed-up).
Usage: experimentalist_k5par.py STEPS MARGIN CELL [--jobs J]   CELL like "11,2:12,2"  (pair (8:3)(9:3) + two bystanders)"""
import sys, time
sys.path.insert(0, "workshop/rounds/081")
from multiprocessing import Pool
from quivermutation import lnaMoves as lm, overlap as ov, nakayama, mutation, pathAlgebra, lines
import theorist_twoby as t

def setup(b1, b2, margin):
    pat = t.PAIR + [b1, b2]
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo); pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + margin + 7
    return L, lm.embedPattern(L, pat, 0), lo, hi

def first_moves(L, rl, lo, hi, margin):
    alg = nakayama.LinearNakayamaAlgebra(L, rl); dual = lm._quiet(pathAlgebra.dualPathAlgebra, alg)
    out = []
    for v in range(max(1, lo - margin), min(L, hi + 1 + margin) + 1):
        for s in (v, -v):
            if lm._quiet(mutation.mutationIsPossibleAtVertex, alg if s > 0 else dual, v): out.append(s)
    return out

def work(args):
    b1, b2, steps, margin, s = args
    L, rl, lo, hi = setup(b1, b2, margin)
    allowed = list(range(max(1, lo - margin), min(L, hi + 1 + margin) + 1))
    start = lines.className(rl); best = {}; seen = {}
    def walk(alg, st, hist):
        if st == 0: return
        dual = lm._quiet(pathAlgebra.dualPathAlgebra, alg)
        for v in allowed:
            for sg in (v, -v):
                if not lm._quiet(mutation.mutationIsPossibleAtVertex, alg if sg > 0 else dual, v): continue
                nx = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [sg])
                if nx is None: continue
                key = lm._stateKey(nx); seq = hist + [sg]
                r = lm.asRelLengths(nx, L)
                if r is not None:
                    nm = lines.className(r)
                    if nm != start and (nm not in best or len(seq) < len(best[nm])): best[nm] = seq
                if seen.get(key, -1) >= st - 1: continue
                seen[key] = st - 1
                walk(nx, st - 1, seq)
    alg = nakayama.LinearNakayamaAlgebra(L, rl)
    if s < 0: pass
    nx = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s])
    if nx is not None:
        r = lm.asRelLengths(nx, L)
        if r is not None:
            nm = lines.className(r)
            if nm != start: best[nm] = [s]
        walk(nx, steps - 1, [s])
    return best

if __name__ == "__main__":
    steps, margin = int(sys.argv[1]), int(sys.argv[2])
    b1, b2 = [tuple(map(int, b.split(","))) for b in sys.argv[3].split(":")]
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
    L, rl, lo, hi = setup(b1, b2, margin)
    fm = first_moves(L, rl, lo, hi, margin); t0 = time.time()
    print("cell", b1, b2, "steps", steps, "margin", margin, "first moves", len(fm), flush=True)
    merged = {}
    with Pool(jobs) as p:
        for i, b in enumerate(p.imap_unordered(work, [(b1, b2, steps, margin, s) for s in fm])):
            for k, v in b.items():
                if k not in merged or len(v) < len(merged[k]): merged[k] = v
            print("done", i + 1, "/", len(fm), "reached so far", len(merged), round(time.time() - t0), flush=True)
    ovs = [ov.maxOverlap([int(c) for c in n]) for n in merged]; base = ov.maxOverlap(rl)
    print("RESULT", b1, b2, "base", base, "reached", len(merged), "min", min(ovs) if ovs else None,
          "lowered", sum(o < base for o in ovs), "secs", round(time.time() - t0, 1), flush=True)
