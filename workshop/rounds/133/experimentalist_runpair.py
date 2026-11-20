"""Round 133 referee items 2-3.  Usage: experimentalist_runpair.py W P
(a) Rules of lnaMoves.VERIFIED_MOVES whose before-pattern has overlap 2 and after-pattern <= 1;
    how many have a gap-1 (3,3) pair in the before-pattern.
(b) The band closure of experimentalist_band.py (doubles + VERIFIED_MOVES), overlap-2 rows tallied
    by (longest overlap>=2 run, (3,3) pair anywhere in the row, (3,3) pair inside a longest run)."""
import sys, time
from collections import defaultdict
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from experimentalist_band import band_rows, edges
from quivermutation import overlap as ov, lnaMoves as lm

def pat_ov(rel):
    return max([s + x - s2 for (s, x), (s2, y) in zip(rel, rel[1:])] + [0])

def has33(rel):
    return any(x == 3 and y == 3 and s2 == s + 1 for (s, x), (s2, y) in zip(rel, rel[1:]))

def runs(rel):                                   # maximal runs linked by overlap >= 2
    out, cur = [], [rel[0]] if rel else []
    for p, q in zip(rel, rel[1:]):
        if p[0] + p[1] - q[0] >= 2: cur.append(q)
        else: out.append(cur); cur = [q]
    if cur: out.append(cur)
    return out

def reach_mode(W, P):                            # `W P --reach`: clean (4,3) forward closure, refined tally
    from experimentalist_band import clean_cells
    L = W + 2 * P; start = [r for r in band_rows(W, P) if clean_cells(r, 4, 3)]
    seen = set(start); fr = list(start)
    while fr:
        for y in edges(fr.pop(), L, W, P, "both"):
            if y not in seen: seen.add(y); fr.append(y)
    tal = defaultdict(int)
    for r in seen:
        if ov.maxOverlap(r) != 2: continue
        rel = lm.relationsOf(list(r)); rs = runs(rel); m = max(map(len, rs))
        tal[(m, has33(rel), any(len(x) == m and has33(x) for x in rs))] += 1
    print(f"reach W {W}: clean (4,3) {len(start)}, reached {len(seen)}; overlap-2 reached rows by"
          " (longest run, (3,3) in row, (3,3) in a longest run):", sorted(tal.items()))

if __name__ == "__main__" and "--reach" in sys.argv:
    reach_mode(int(sys.argv[1]), int(sys.argv[2])); sys.exit()

if __name__ == "__main__":
    lowr = [d for d in lm.VERIFIED_MOVES if pat_ov(list(d[1])) == 2 and pat_ov(list(d[2])) <= 1]
    print("rules", len(lm.VERIFIED_MOVES), "lowering own overlap 2 -> <=1:", len(lowr),
          "with gap-1 (3,3) pair before:", sum(has33(list(d[1])) for d in lowr))
    W, P = int(sys.argv[1]), int(sys.argv[2]); L = W + 2 * P; t0 = time.time()
    rows = band_rows(W, P); back = defaultdict(list)
    for r in rows:
        for r2 in edges(r, L, W, P, "both"):
            back[r2].append(r)
    seen = {r for r in rows if ov.maxOverlap(r) <= 1}; fr = list(seen)
    while fr:
        for y in back[fr.pop()]:
            if y not in seen: seen.add(y); fr.append(y)
    tal = defaultdict(lambda: [0, 0])
    for r in rows:
        if ov.maxOverlap(r) != 2: continue
        rel = lm.relationsOf(list(r)); rs = runs(rel); m = max(map(len, rs))
        k = (m, has33(rel), any(len(x) == m and has33(x) for x in rs))
        tal[k][0] += 1; tal[k][1] += r in seen
    print(f"W {W} P {P} secs {time.time()-t0:.1f}")
    print("(longest run, (3,3) in row, (3,3) in a longest run) | rows | reach <= 1")
    for k, (n, l) in sorted(tal.items()):
        print(k, "|", n, "|", l)
