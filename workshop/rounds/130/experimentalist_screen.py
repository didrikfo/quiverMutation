"""T7/H-010 round 130: screen for an interior start-gap-2 control.
Usage: experimentalist_screen.py DEPTH [--pairs "8,4:10,3;..."] [--bys "2..16" ] [--mmax 6] [--two] [--clean]
--two: two bystanders instead of one.  --clean: keep only cells whose base equals the pair's own
overlap and no bystander starts at a pair start +-1 (so the cell is not a gap-1 run in disguise).
Each cell = a start-gap-2 pair (default E-330's four interior shapes) plus ONE bystander (s, m),
admissible, placed far from both ends (L = hi + 9, min start >= 7).  Base = maxOverlap of the cell.
Cheap pre-screen: closure under doubleMutation.rewritesOf (each rewrite = 2 mutations), DEPTH
rewrites, keeping only rewrites whose mutation vertices lie in the margin-1 window
[lo - 1, hi + 2] that localMutationSequences(margin=1) would allow.  Prints every cell that
reaches maxOverlap < base, with its sequence; and a count line.  A screen, not a verdict:
cells it misses may still lower by single mutations (engine check: rounds/127 experimentalist_gap2fix.py)."""
import sys
from quivermutation import lnaMoves as lm, overlap as ov, doubleMutation as dm

def cells(pairs, bys, mmax, two=False):
    one = [(s, m) for s in bys for m in range(2, mmax + 1)]
    for pair in pairs:
        starts = {p for p, _ in pair}
        for i, b in enumerate(one):
            for c in (one[i + 1:] if two else [None]):
                extra = [b] + ([c] if c else [])
                if len({s for s, _ in extra} | starts) != len(extra) + len(starts):
                    continue
                yield pair, tuple(extra), sorted(list(pair) + extra)

def screen(pat, depth):
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo); pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + 9
    rl = lm.embedPattern(L, pat, 0)
    if rl is None:
        return None
    win = set(range(lo - 1, hi + 3))
    base = ov.maxOverlap(rl)
    seen = {tuple(rl): []}; frontier = [tuple(rl)]
    for _ in range(depth):
        nxt = []
        for row in frontier:
            for row2, seq in dm.rewritesOf(L, row, allowSource=False):
                if row2 is None or not all(abs(v) in win for v in seq) or row2 in seen:
                    continue
                seen[row2] = seen[row] + list(seq); nxt.append(row2)
        frontier = nxt
    best = min(seen, key=lambda r: ov.maxOverlap(list(r)))
    return dict(L=L, pat=pat, base=base, reached=len(seen) - 1,
                min=ov.maxOverlap(list(best)), best="".join(map(str, best)), seq=seen[best])

if __name__ == "__main__":
    depth = int(sys.argv[1])
    a = sys.argv
    pairs = ([tuple(tuple(map(int, r.split(","))) for r in c.split(":")) for c in a[a.index("--pairs") + 1].split(";")]
             if "--pairs" in a else [((8, 4), (10, 3)), ((8, 5), (10, 4)), ((8, 4), (10, 4)), ((8, 5), (10, 5))])
    b0, b1 = map(int, a[a.index("--bys") + 1].split("..")) if "--bys" in a else (2, 16)
    mmax = int(a[a.index("--mmax") + 1]) if "--mmax" in a else 6
    two, clean = "--two" in a, "--clean" in a
    n = adm = low = 0
    for pair, by, pat in cells(pairs, range(b0, b1 + 1), mmax, two):
        n += 1
        r = screen(pat, depth)
        if r is None:
            continue
        pov = pair[0][0] + pair[0][1] - pair[1][0]
        if clean and (r["base"] != pov or any(abs(s - p) == 1 for s, _ in by for p, _ in pair)):
            continue
        adm += 1
        if r["min"] < r["base"]:
            low += 1
            print("LOWERS", pair, "by", by, "pairOv", pov, r, flush=True)
    print("cells", n, "admissible" + (" and clean" if clean else ""), adm, "lowering", low)
