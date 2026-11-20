"""T7/H-010 round 133: exhaustive interior closure on a band, all bystander configurations at once.
Usage: experimentalist_band.py W P [--moves dbl|rules|both]   (default both)
Every LNA whose relations lie in the band of W vertices [P+1, P+W] of A_L, L = W + 2P (so P free
vertices at each end), is a node.  Edges: the library's derived moves restricted to the band --
`doubleMutation.rewritesOf` (allowSource False; 'dbl') and/or the floating verified table
`lnaMoves.VERIFIED_MOVES` ('rules'; no anchored/end rules, no edge moves, no free move) -- kept
only when the result is again in the band.  Backward BFS from the rows of maxOverlap <= m gives,
for each m, every band row that some interior move sequence (any length) takes to overlap <= m.
Reported per start-gap-2 shape (a,b): clean cells = rows containing consecutive relations
(s,a),(s+2,b) with maxOverlap == overlap of the pair, no relation starting at s-1, s+1, s+3;
how many of them reach maxOverlap < base.  Positive controls: (5,4), (6,5) at overlap 3, 4
(E-336 found interior lowering for (8,5)(10,4)), and the run 333 (overlap 2, gap 1).
A statement about this move set only; single-mutation routes are the engine's job."""
import sys, time
from collections import defaultdict
from quivermutation import nakayama, overlap as ov, doubleMutation as dm, lnaMoves as lm

def band_rows(W, P):
    return [(0,) * P + tuple(r) + (0,) * P for r in nakayama.allRelationLengths(W)]

def in_band(row, W, P):
    return all(v == 0 for v in row[:P]) and all(
        v == 0 or i + 1 + v <= P + W for i, v in enumerate(row))

def edges(row, L, W, P, moves):
    out = set()
    if moves in ("dbl", "both"):
        out |= {r for r, _ in dm.rewritesOf(L, row, allowSource=False)}
    if moves in ("rules", "both"):
        out |= set(lm.rewritesOf(L, list(row), lm.VERIFIED_MOVES))
    return [r for r in out if r != row and in_band(r, W, P)]

def clean_cells(row, a, b):
    rel = lm.relationsOf(list(row)); starts = {s for s, _ in rel}
    for (s, x), (s2, y) in zip(rel, rel[1:]):
        if x == a and y == b and s2 == s + 2 and ov.maxOverlap(row) == a - 2 \
                and not starts & {s - 1, s + 1, s + 3}:
            return True
    return False

if __name__ == "__main__":
    W, P = int(sys.argv[1]), int(sys.argv[2]); L = W + 2 * P
    moves = sys.argv[sys.argv.index("--moves") + 1] if "--moves" in sys.argv else "both"
    t0 = time.time()
    rows = band_rows(W, P)
    back = defaultdict(list); nedge = 0
    for r in rows:
        for r2 in edges(r, L, W, P, moves):
            back[r2].append(r); nedge += 1
    mo = {r: ov.maxOverlap(r) for r in rows}
    reach = {}                                   # m -> set of rows that reach maxOverlap <= m
    for m in (1, 2, 3):
        seen = {r for r in rows if mo[r] <= m}; fr = list(seen)
        while fr:
            x = fr.pop()
            for y in back[x]:
                if y not in seen:
                    seen.add(y); fr.append(y)
        reach[m] = seen
    print(f"W {W} P {P} L {L} moves {moves} rows {len(rows)} edges {nedge} secs {time.time()-t0:.1f}")
    print("shape | pair overlap | clean cells | isolated(moves 0) | reach overlap < base | example")
    for a, b in [(4, 3), (4, 4), (5, 4), (5, 5), (6, 5)]:
        cells = [r for r in rows if clean_cells(r, a, b)]
        if not cells:
            continue
        m = a - 3
        low = [r for r in cells if r in reach[m]]
        iso = sum(1 for r in cells if not edges(r, L, W, P, moves))
        ex = "".join(map(str, low[0])) if low else "-"
        print(f"({a},{b}) | {a-2} | {len(cells)} | {iso} | {len(low)} | {ex}")
    ov2 = [r for r in rows if mo[r] == 2]
    print("all overlap-2 rows", len(ov2), "reach overlap <= 1:", sum(r in reach[1] for r in ov2))
    # which overlap-2 rows reach almost separate: tally by the multiset of relation lengths >= 3
    tal = defaultdict(int)
    for r in ov2:
        if r in reach[1]:
            tal[tuple(sorted(v for v in r if v >= 3))] += 1
    print("overlap-2 rows reaching <=1, by lengths>=3 multiset (top 12):",
          sorted(tal.items(), key=lambda kv: -kv[1])[:12])
    has4 = [r for r in ov2 if any(v >= 4 for v in r)]
    print("overlap-2 rows with a relation of >= 4 arrows:", len(has4),
          "reach <= 1:", sum(r in reach[1] for r in has4))
    # by the set of overlap-2 consecutive pair types (a, b): how many rows, how many reach <= 1
    def pairtypes(r):
        rel = lm.relationsOf(list(r))
        return tuple(sorted({(x, y) for (s, x), (s2, y) in zip(rel, rel[1:])
                             if s + x - s2 == 2}))
    by = defaultdict(lambda: [0, 0])
    for r in ov2:
        k = pairtypes(r); by[k][0] += 1; by[k][1] += r in reach[1]
    print("overlap-2 pair types | rows | reach <= 1")
    for k, (n, l) in sorted(by.items(), key=lambda kv: -kv[1][0]):
        if l or n >= 20:
            print(k, "|", n, "|", l)
    no33 = [r for r in ov2 if (3, 3) not in pairtypes(r)]
    print("overlap-2 rows with no (3,3) gap-1 pair:", len(no33), "reach <= 1:",
          sum(r in reach[1] for r in no33))
    # by the longest run of relations linked by overlaps >= 2, and whether that run holds a (3,3) pair
    def longest_run(r):
        rel = lm.relationsOf(list(r)); best = cur = 1
        for (s, x), (s2, y) in zip(rel, rel[1:]):
            cur = cur + 1 if s + x - s2 >= 2 else 1; best = max(best, cur)
        return best
    rr = defaultdict(lambda: [0, 0])
    for r in ov2:
        k = (longest_run(r), (3, 3) in pairtypes(r)); rr[k][0] += 1; rr[k][1] += r in reach[1]
    print("(longest overlap>=2 run, has (3,3)) | rows | reach <= 1")
    for k, (n, l) in sorted(rr.items()):
        print(k, "|", n, "|", l)
