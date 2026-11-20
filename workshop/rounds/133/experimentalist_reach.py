"""T7/H-010 round 133: what the clean (4,3) cells of a band reach under the library's derived moves.
Usage: experimentalist_reach.py W P
Same band, move set and cleanness as experimentalist_band.py (double + VERIFIED_MOVES, any length, rows
kept in the band).  Forward closure from all clean (4,3) cells at once; tallies the reached rows by
(maxOverlap, longest overlap>=2 run, has a gap-1 (3,3) pair) -- the category the band run found to
decide lowering at overlap 2 (only longest run 3 with a (3,3) pair ever reaches overlap <= 1)."""
import sys, time
from collections import Counter
from quivermutation import overlap as ov, lnaMoves as lm
from experimentalist_band import band_rows, edges, clean_cells

def cat(r):
    rel = lm.relationsOf(list(r)); best = cur = 1; has33 = False
    for (s, x), (s2, y) in zip(rel, rel[1:]):
        o = s + x - s2
        cur = cur + 1 if o >= 2 else 1; best = max(best, cur)
        has33 |= (x, y, o) == (3, 3, 2)
    return (ov.maxOverlap(r), best, has33)

if __name__ == "__main__":
    W, P = int(sys.argv[1]), int(sys.argv[2]); L = W + 2 * P; t0 = time.time()
    start = [r for r in band_rows(W, P) if clean_cells(r, 4, 3)]
    seen = set(start); fr = list(start)
    while fr:
        x = fr.pop()
        for y in edges(x, L, W, P, "both"):
            if y not in seen:
                seen.add(y); fr.append(y)
    print(f"W {W} P {P} clean (4,3) cells {len(start)} reached rows {len(seen)} secs {time.time()-t0:.1f}")
    print("(maxOverlap, longest run, has (3,3)) | rows reached")
    for k, n in sorted(Counter(cat(r) for r in seen).items()):
        print(k, "|", n)
