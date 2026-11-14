"""T7/H-010 round 127: gap-2 pair cells after the E-313 ValueError fix (lnaMoves._copy keeps arrowRels).
Usage: experimentalist_gap2fix.py STEPS MARGIN CAPSECS "s,l:s,l;..." [--end] [--jobs J]
Same harness as rounds/107 experimentalist_chain.run (pattern shifted to min start 7, L = hi+margin+7,
mutations only within MARGIN of the pattern), except --end plants the pattern at start 1 so the
left end of the quiver is inside the mutation window (H-011's end mechanism: the planted control).
Per-cell wall cap; a cap or an exception is printed as such, never as a verdict.
Prints: pattern, base overlap, #LNAs reached, min overlap, #lowered, a lowest LNA and its sequence, secs."""
import sys, signal, time
from multiprocessing import Pool
from quivermutation import lnaMoves as lm, overlap as ov

def _alarm(*a): raise TimeoutError

def run(a):
    pat, steps, margin, cap, end = a
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(cap)
    t0 = time.time()
    try:
        pat = list(pat)
        lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
        off = (1 - lo) if end else max(0, 7 - lo)
        pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
        L = hi + margin + 7
        rl = lm.embedPattern(L, pat, 0)
        if rl is None: return (a[0], "inadmissible")
        reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
        base = ov.maxOverlap(rl)
        ovs = {n: ov.maxOverlap([int(c) for c in n]) for n in reached}
        best = min(ovs, key=ovs.get) if ovs else None
        return (a[0], "end" if end else "interior", "L", L, "base", base, "reached", len(reached),
                "min", ovs[best] if ovs else None, "lowered", sum(o < base for o in ovs.values()),
                "best", best, reached.get(best), "secs", round(time.time() - t0, 1))
    except TimeoutError:
        return (a[0], "CAP", cap, "s: undecided (not a verdict)")
    except Exception as e:
        return (a[0], "ERROR", type(e).__name__, str(e)[:90])
    finally:
        signal.alarm(0)

if __name__ == "__main__":
    steps, margin, cap = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    end = "--end" in sys.argv
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
    cells = [tuple(tuple(map(int, r.split(","))) for r in c.split(":")) for c in sys.argv[4].split(";")]
    print("steps", steps, "margin", margin, "cap", cap, "end" if end else "interior", flush=True)
    with Pool(jobs) as p:
        for r in p.imap_unordered(run, [(c, steps, margin, cap, end) for c in cells]):
            print(*r, flush=True)
