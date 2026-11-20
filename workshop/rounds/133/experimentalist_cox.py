"""T7/H-010 round 133: a derived-invariant obstruction for overlap-2 cells.
Usage: experimentalist_cox.py count LMAX           -- number of almost separate rows per length
       experimentalist_cox.py band W P             -- band rows of experimentalist_band.py, length L = W + 2P
       experimentalist_cox.py cells "c;c;..."      -- engine cells (gap2fix strings), embedded as the engine does
An LNA is almost separate iff maxOverlap <= 1; every such LNA is derived equivalent to a quipu (Ladkani,
GLOSSARY "almost separate").  S_L = the Coxeter polynomials of the almost separate LNAs of length L.
The Coxeter polynomial is a derived invariant, so an LNA whose polynomial is not in S_L is derived
equivalent to NO almost separate LNA: no sequence of derived equivalences (tilting mutations, any number,
interior or not) takes it to overlap <= 1.  A polynomial IN S_L settles nothing (F-010).
band: tallies clean (4,3), (4,4), (5,4) cells (cleanness as in experimentalist_band.py) and the band rows the
library moves take to overlap <= 1 (the YES rows: their polynomial must be in S_L, a consistency control)."""
import sys, time
from quivermutation import coxeterTables as ct, overlap as ov, lnaMoves as lm

def almost_separate(L):
    """All relation-length rows of length L with every consecutive overlap <= 1 (relations of >= 2 arrows)."""
    out = []
    def rec(row, minstart, lastend):          # lastend = end vertex of the last relation (0 if none)
        out.append(tuple(row))
        for s in range(minstart, L):
            for k in range(2, L - s + 1):     # relation s -> s + k, vertices 1..L
                e = s + k
                if lastend and lastend - s > 1:
                    continue                  # overlap with previous relation > 1
                if e <= lastend:
                    continue                  # ends must increase
                row[s - 1] = k
                rec(row, s + 1, e)
                row[s - 1] = 0
    rec([0] * L, 1, 0)
    return out

def keyset(L):
    t = time.time(); rows = almost_separate(L)
    bad = sum(1 for r in rows if ov.maxOverlap(r) > 1)
    S = {ct.lnaCoxeterKey(L, r) for r in rows}
    print(f"L {L} almost separate rows {len(rows)} (maxOverlap > 1 among them: {bad}) polynomials {len(S)} "
          f"secs {time.time()-t:.1f}", flush=True)
    return S

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "count":
        for L in range(4, int(sys.argv[2]) + 1):
            print(L, len(almost_separate(L)), flush=True)
    elif mode == "band":
        sys.path.insert(0, __file__.rsplit("/", 1)[0])
        from experimentalist_band import band_rows, edges, clean_cells
        from collections import defaultdict
        W, P = int(sys.argv[2]), int(sys.argv[3]); L = W + 2 * P
        S = keyset(L); rows = band_rows(W, P); t = time.time()
        for a, b in [(4, 3), (4, 4), (5, 4)]:
            cells = [r for r in rows if clean_cells(r, a, b)]
            out = sum(1 for r in cells if ct.lnaCoxeterKey(L, r) not in S)
            print(f"clean ({a},{b}) cells {len(cells)} polynomial not in S_L (provably never overlap <= 1): {out}",
                  flush=True)
        # YES control: band rows of overlap 2 that the library moves take to overlap <= 1
        back = defaultdict(list)
        for r in rows:
            for r2 in edges(r, L, W, P, "both"):
                back[r2].append(r)
        seen = {r for r in rows if ov.maxOverlap(r) <= 1}; fr = list(seen)
        while fr:
            x = fr.pop()
            for y in back[x]:
                if y not in seen:
                    seen.add(y); fr.append(y)
        yes = [r for r in seen if ov.maxOverlap(r) == 2]
        print(f"YES control: overlap-2 band rows reaching <= 1: {len(yes)}, polynomial not in S_L: "
              f"{sum(1 for r in yes if ct.lnaCoxeterKey(L, r) not in S)}")
        o2 = [r for r in rows if ov.maxOverlap(r) == 2]
        print(f"all overlap-2 band rows {len(o2)}, polynomial not in S_L: "
              f"{sum(1 for r in o2 if ct.lnaCoxeterKey(L, r) not in S)}  secs {time.time()-t:.1f}")
    elif mode == "cells":
        cells = [[tuple(map(int, x.split(","))) for x in c.split(":")] for c in sys.argv[2].split(";")]
        cache = {}
        for pat in cells:
            lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
            off = max(0, 7 - lo); pat2 = [(s + off, k) for s, k in pat]; L = hi + off + 1 + 7
            row = tuple(lm.embedPattern(L, pat2, 0))
            if L not in cache:
                cache[L] = keyset(L)
            print(pat, "L", L, "overlap", ov.maxOverlap(row), "in S_L" if ct.lnaCoxeterKey(L, row) in cache[L]
                  else "NOT in S_L (no derived equivalence reaches overlap <= 1)", flush=True)
