"""k(c), d for prefix+x words at one n: classes of held offsets (orbit, no mirror join), pair sums s, k = n - s, d = singleton count.
  timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_kd.py N PREFIX LO-HI [LIMIT]
Output one line per word: n word hi classes | sums | k | d | status  (CAP = orbit walk hit the limit: undecided, not a verdict)."""
import sys, time
import batch
from quivermutation import freeMoves
n, pre = int(sys.argv[1]), sys.argv[2]
lo, hi_x = map(int, sys.argv[3].split("-"))
limit = int(sys.argv[4]) if len(sys.argv) > 4 else 300000
R = freeMoves.REDUCED
for x in range(lo, hi_x + 1):
    word = "%s%d" % (pre, x)
    t0 = time.time()
    offs = [o for o in range(n) if batch._rowFor(n, word, o) is not None]
    if not offs:
        print(n, word, "no placement", flush=True); continue
    done, classes, cap, sizes = set(), [], False, []
    for o in offs:
        if o in done: continue
        w = freeMoves.orbitReport(n, batch._rowFor(n, word, o), free=R, limit=limit)
        if not w.closed: cap = True
        held = [p for p in offs if freeMoves._startOf(batch._rowFor(n, word, p), R) in w.rows]
        done.update(held); classes.append(held); sizes.append(len(w.rows))
    pairs = [c for c in classes if len(c) == 2]
    sums = sorted(set(c[0] + c[1] for c in pairs))
    multi = [c for c in classes if len(c) > 2]
    k = (n - sums[0]) if len(sums) == 1 else None
    d = sum(1 for c in classes if len(c) == 1)
    print(n, word, "hi", offs[-1], "classes", classes, "sizes", sizes, "| sums", sums, "| k", k, "| d", d,
          "| big-classes", len(multi), "| CAP" if cap else "| closed", "%.0fs" % (time.time() - t0), flush=True)
