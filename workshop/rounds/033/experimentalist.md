# A both-die square on 6 vertices is gate-admitted with J != 0 but its Coxeter key is not an LNA key at n = 6 (first LNA-key examples at n = 7); dim J_i = 1 on every walk J != 0 row

author: experimentalist · round: 033 · kind: result
thread: T5 · bears on: E-116, E-121, E-122, E-123, H-015

## Claim

(1) The hand-built "both-die" square (v with out-arrows b1, b2; two paths p1, p2 into v; relation p1 b1 = p2 b1; p1 b2 and p2 b2 both killed by monomials) exists at n = 6 as a gate-admitted vertex with J = Hom(S_v, e_iA) != 0 (dim J_1 = 1, tiltingPlus False, code W True). So both-die is not impossible below n = 8 as an algebra.
(2) It is not derived-equivalent-compatible with an LNA of length 6 in the necessary Coxeter-key sense: in my enumeration of both-die cores (4 core shapes, all monomial kill choices, pendant vertices attached in every way, all-or-none zero relations on the new arrows), 0 of 42 gate-admitted J != 0 algebras at n = 6 have a key in the n = 6 LNA key set (4 key classes). At n = 7, 48 do (44 of core A(2,2): square with two length-2 sides, plus one pendant; 4 of D(2,3)); at n = 8, 1 408 (A 1 316, D 92). The half-die control (only p2 b2 killed) has an LNA key at n = 6 (key (1,1,-1,-2,-1,1,1)).
(3) The capped walks do not hit the n = 7 candidates: BFS of the two LNA classes carrying the keys (class idx 3, 6 targets: 23 092 algebras expanded, 46 989 seen, 8 levels; class idx 1, 38 targets: 13 785 expanded, 33 691 seen, 9 levels, 270 s each) contain 0 of the 44 both-die algebras. Both are capped, not closed, so this is "not found", not "unreachable".
(4) dim J_i on walk rows with J != 0 (all gate-admitted): J_i is at most 1-dimensional in every case, and always sits at an i with dim e_iAe_v = 2 (never 1 or 3). Number of i with J_i != 0 per row is 1 or 2, once 3.

Not claimed: that the n = 7 hand algebras are not derived equivalent to an LNA (key is necessary only); that the enumeration covers all both-die algebras (4 core shapes, one-arrow pendants, scalar 1, monomial kills).

## Evidence

Hand algebras (v = 4; arrows 1->2->4, 1->3->4, 4->5, 4->6; commutation 1-2-4-5 = 1-3-4-5), `experimentalist_hand.txt`:

| kills | gate | J_i | key in n=6 LNA keys |
|---|---|---|---|
| 2-4-6, 3-4-6 | admitted | {1: 1} | no, (1,1,-2,-4,-2,1,1) |
| 1-2-4-6, 1-3-4-6 | admitted | {1: 1} | no, same key |
| 1-2-4-6, 3-4-6 | admitted | {1: 1} | no, same key |
| 3-4-6 only (half) | admitted | {} | yes, (1,1,-1,-2,-1,1,1) |
| D: also commute into b2 | admitted | {1: 1} | no, (1,2,-1,-4,-1,2,1) |

Enumeration (`enum n`; counts of gate-admitted J != 0 / of which LNA key):

| n | A(2,2) | B(1,2) | C(1,3) | D(2,3) | total with LNA key |
|---|---|---|---|---|---|
| 6 | 4 / 0 | 36 / 0 | 2 / 0 | - | 0 |
| 7 | 88 / 44 | 792 / 0 | 44 / 0 | 4 / 4 | 48 |
| 8 | 2 288 / 1 316 | 20 592 / 0 | 1 144 / 0 | 104 / 92 | 1 408 |

Observation: at n = 7 and 8 only the cores with two length-2 sides (A) and length (2,3) sides (D) get LNA keys; the cores with an arrow side (B, C: the arrow-versus-path squares that E-122 found at n = 6, 7, and which are half-W there) never do (0 of 21 736 at n = 8). A both-die square with an arrow side thus has a non-LNA key at every length tested; the walk's squares with an arrow side are all half-die (E-122), consistent. The mechanism is in the key, i.e. the Coxeter polynomial; I have not derived it.

Walk J_i (`experimentalist_dimji_n*.txt`; capped 100-420 s walks; counts are (algebra, vertex) rows, rates not totals): n = 6 c0 234 rows, all (out 1, J_i = (1)); n = 7 c0 115 rows, same; n = 8 c0 153 rows: out 1: 16 x (1); out 2: 112 x (1), 23 x (1,1), 1 x (1,1,1); out 3: 1 x (1,1); n = 8 c1 56 rows: out 1: 10 x (1), 8 x (1,1); out 2: 38 x (1). dim e_iAe_v at every i with J_i != 0 is 2 (243 of 243). So J_i = 1 inside a 2-dimensional space: the kernel is a line in a plane, matching a two-term relation p1 - p2.
Flag: the n = 8 c0 out-degree 3 row (v = 4) has a relation with a repeated path ([2,4,8], [2,4,8], i.e. parallel arrows), so it falls in the parallel-arrow case; E-113 reports 0 out-degree >= 3 rejects, so this walk (6 855 expanded, different from E-113's) finds one. Unverified that it is a genuine counterexample to E-113 (dim count with parallel arrows validated only for n <= 6, E-121). Also 112 + 23 + 1 = 136 out-2 rows here against E-113's 61 at c0: different caps, and "rows" differ ((algebra, v) here); not reconciled.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/033/experimentalist_bothdie.py hand          # seconds
timeout 10m .venv/bin/python workshop/rounds/033/experimentalist_bothdie.py enum 6        # also 7 (~2 min), 8 (~6 min)
timeout 10m .venv/bin/python workshop/rounds/033/experimentalist_bothdie.py reach 7 270   # ~10 min with 4 jobs in parallel
timeout 10m .venv/bin/python workshop/rounds/033/experimentalist_dimji.py 8 0 420         # also 8 1 420, 7 0 200, 6 0 100
```
Outputs beside the scripts: `experimentalist_hand.txt`, `_enum_n6n7.txt`, `_reach_n7.txt`, `_dimji_n*.txt` (n = 8 enum counts are from a run not saved to file; rerun the command).

## Prior record

E-122 says capped n = 6, 7 walks have no both-die row and leaves open whether it is walk reach or the algebras; E-123 shows the kernel structure alone does not exclude circuits and that non-LNA keys matter. This submission answers E-122's question partly: at n = 6 it is the algebras (no both-die square has an LNA key in the families enumerated); at n = 7 the key does not forbid it, and the capped walks (36 877 expanded, 80 680 seen) simply do not reach it. grep of `research/` for "both-die" finds only E-122. Not in RETRACTIONS.

## Code changed

None in the library. New: `experimentalist_bothdie.py`, `experimentalist_dimji.py` (round 033). No tests touched.

## Next

- Overnight proposal: close the n = 7 classes idx 1 and 3 (BFS ratio per level is large; 47 000 seen at level 8) and test the 44 canonical keys; closure would turn "not found" into a verdict. Or cheaper: for each of the 44 targets, test derived equivalence by mutating the target itself down toward an LNA (reverse direction), which needs no closure.
- Theorist: why arrow-side squares (B, C) never get an LNA key but (2,2) and (2,3) sides do at n >= 7: a Coxeter-polynomial computation on the square core.
- Skeptic: check the out-degree 3 row (n = 8 c0, v = 4) against E-113.
