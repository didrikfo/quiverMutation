"""T7/H-010 round 133: pick clean two-bystander (8,4)(10,3) cells for the engine, and find the shortest
library-move route to overlap <= 1 for a given cell (used for the non-clean 4033 control).
Usage: experimentalist_pick.py pick [N] [SEED]   -- N clean, non-vacuous (screen depth 3) cells, bystander
          starts 5..15, lengths 2..5, both bystanders meeting the margin-1 window [7, 13] (so the engine
          can touch them); random sample with SEED; prints the gap2fix cell string.
       experimentalist_pick.py route "8,4:10,3:11,3" [DEPTH]  -- forward BFS from the E-330-embedded cell
          (min start >= 7, L = hi + 8) under doubleMutation (allowSource False) + VERIFIED_MOVES, mutation
          vertices kept within margin 1 of the pattern span; prints the shortest sequence to overlap < base."""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "130"))
from experimentalist_screen import screen, cells
from quivermutation import lnaMoves as lm, overlap as ov, doubleMutation as dm

def pick(n, seed):
    pair = ((8, 4), (10, 3))
    out = []
    for p, by, pat in cells([pair], range(5, 16), 5, two=True):
        if any(abs(s - q) == 1 for s, _ in by for q, _ in pair):
            continue
        if not all(s <= 14 and s + m >= 7 for s, m in by):     # bystander inside/at the window
            continue
        r = screen(pat, 3)
        if r is None or r["base"] != 2 or r["reached"] == 0:
            continue
        out.append(pat)
    random.Random(seed).shuffle(out)
    print("eligible", len(out))
    sel = sorted(out[:n])
    print(";".join(":".join(f"{s},{k}" for s, k in c) for c in sel))

def route(cell, depth):
    pat = [tuple(map(int, r.split(","))) for r in cell.split(":")]
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo); pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + 1 + 7
    rl = tuple(lm.embedPattern(L, pat, 0)); base = ov.maxOverlap(rl)
    win = set(range(lo - 1, hi + 3))
    seen = {rl: []}; fr = [rl]
    for d in range(depth):
        nxt = []
        for row in fr:
            moves = [(r, s) for r, s in dm.rewritesOf(L, row, allowSource=False)]
            for desc, ws in lm.matchingWindows(L, list(row), lm.VERIFIED_MOVES):
                a = lm.applyAt(L, list(row), desc, ws)
                if a is not None:
                    moves.append((tuple(a[0]), list(a[1])))
            for r, s in moves:
                if r in seen or not all(abs(v) in win for v in s):
                    continue
                seen[r] = seen[row] + list(s); nxt.append(r)
                if ov.maxOverlap(r) < base:
                    print("L", L, "base", base, "->", "".join(map(str, r)), "overlap", ov.maxOverlap(r),
                          "mutations", len(seen[r]), seen[r]); return
        fr = nxt
    print("L", L, "base", base, "no route to overlap <", base, "within", depth, "moves;", len(seen), "rows")

if __name__ == "__main__":
    if sys.argv[1] == "pick":
        pick(int(sys.argv[2]) if len(sys.argv) > 2 else 16, int(sys.argv[3]) if len(sys.argv) > 3 else 133)
    else:
        route(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4)
