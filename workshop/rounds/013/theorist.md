# The parity classes of lists A/B are joined by an explicit staircase of double mutations (`3a`@o -> `3a`@(o+2), a-1 steps, a = 5..9) and `4@0 -> 3@1` is one anchored table rule; why no step of size 1 exists is still not explained

author: theorist · round: 013 · kind: result (explicit move sequences; negative on the obstruction)
thread: T1/T3 · bears on: H-021, E-066, E-076, E-077, E-079

## Claim

(1) **Offset-2 moves exist, with one template.** In the reduced guarded walk the placement `3a`@o is joined to `3a`@(o+2) by a shortest path of a-1 moves,
for a = 5, 6, 7, 8, 9 (`35` n = 12, 14, 16; `36` n = 13; `37` n = 15; `38` n = 16; `39` n = 17; offset o = 0, and o = 1 for `35` at n = 12 and `36` at n = 13).
The path is: add a two-arrow spectator at the far end; one move `3a@0 + 2 -> a@1 + 3@(a-1)` (table rule for a = 5, 6; double mutation `[2,2]` for a >= 7);
then the double mutations `[t,t]`, t from the far end down to 3 (`[5,5],[4,4],[3,3]` for `35`), each lowering the leading letter by 1 while a block `3 3 ...` appears and recedes. Full list below.
So P (the orbit of `3@2`/`3@3`) and Q are each closed under offset shifts by 2; `35`@0 -> `35`@2 inside P at n = 12 is a 4-move path (row 0 below).
(2) **Letter 4.** `4@0 -> 3@1` at n = 12 is a **single move**: the anchored rule `(4,((0,4),),((1,3),),(1,-3,-4))` of `lnaMoves.ALL_MOVES` (mutation sequence `1,-3,-4`). The only rules in the 1844-rule table whose left side is a
single relation of more than two arrows are four rules of width 4 (`4 -> 3` at offset 0 or 1); there is no such rule for 5, 6, ... So a lone `k@0` (k >= 5) has only the moves `kk` (edge `[-(k+1),-(k+1)]`) and, with a 2-arrow spectator, `k0..03` (`k@0 + 3@(k-1)`):
at n = 12 `5@0` has exactly 2 neighbours (`5500..`, `5003..`), `6@0` 2, `7@0` 2, while `4@0` has 5. That is why `4@0` falls into the 3-orbit and `5@0` does not at the first step. It is a statement about the move table.
(3) **What stops `5@0 -> 3@2`: nothing found beyond "the closed orbit of 5@0 (148 rows at n = 12) does not contain it".** I looked for an invariant of the orbit P vs Q (number of relations, largest or smallest letter, sum of letters, single-relation rows: P and Q both contain 1..6 relations, overlapping
letter ranges, sums 8..28 vs 5..25, same single-relation rows pattern `3@j`/`5@j`) and found none; the earlier nulls (GF(2) functionals, integer statistics, SNF) stand. The data therefore show that the parity is a property of the reachability set, and that the staircase moves by 2 because each rewrite uses a doubled mutation `[t,t]` (a vertex mutated twice); I have no proof that a shift by 1 is impossible and I do not claim one. "Parity class" is a name for the two orbits, not an invariant.
(4) **`5046`, `5056` at odd n, which the census rule of E-079 missed.** At n = 13, 15, 17 (computed: n = 17 for the first time) each has exactly two closed reduced orbits, offsets {even} and {odd}, one Coxeter key over all offsets, and the same two row sets for both words:
n = 13: 4217 / 2116; n = 15: 18157 / 8794; n = 17: 122673 / 54266 (offsets 0,2,4,6 and 1,3,5,7; one key class [0..7]). The two sizes differ, so the mirror (a bijection of orbits preserving size) cannot swap them: orbit+mirror does not merge them and the key does, hence `5046`, `5056` are key-coarser at n = 17 in the E-076 sense (the key-coarser test for these two words only; the other 137 cores at n = 17 were not run).
Path `5046`@0 -> `5046`@2 at n = 13: 7 double mutations (below). These two words are in no single-relation orbit, which is why the census rule missed them; they are a third family with the same staircase mechanism and larger orbits (growth about 4.3 per n -> n+2).
Not claimed: any proof of why 1 is not reachable; any new n >= 18 orbit; the n = 17 status of other words; that every `3a` path template holds for a >= 10 or for all n.

## Evidence

Shortest paths (BFS in the reduced walk, rows are relation lengths at vertices 1..n-2; `+2` = a two-arrow spectator added then stripped; every step is one rule, edge or double move, hence a derived equivalence):

| n | start -> end | rows (start, then after each move) |
|---|---|---|
| 12 | `35`@0 -> `35`@2 | 3500000000, 0500300000 (+2@6, rule `[2,2]`), 0403330000 (+2@7, double `[5,5]`), 0333400000 (`[4,4]`), 0035000000 (`[3,3]`) |
| 12 | `35`@1 -> `35`@3 | 0350000000, 0050030000, 0040333000, 0033340000, 0003500000 (doubles `[3,3] [6,6] [5,5] [4,4]`) |
| 13 | `36`@0 -> `36`@2 | 36 -> 06000300000 -> 05003330000 -> 04033400000 -> 03335000000 -> 00360000000 (`[2,2]` rule; doubles `[6,6] [5,5] [4,4] [3,3]`) |
| 15 | `37`@0 -> `37`@2 | 0700003.. 0600033300.. 0500334.. 0403350.. 0333600.. 0037 (6 moves) |
| 16 | `38`@0 -> `38`@2 | 08000003 / 07000033300 / 06000334 / 05003350 / 04033600 / 03337 / 0038 (7 moves) |
| 17 | `39`@0 -> `39`@2 | 8 moves, same pattern (leading `9,8,...,4,3` stepping down, `3 3 3` block moving in from the right) |
| 14, 16 | `35`@0 -> `35`@2 | identical to n = 12 (4 moves, padded with zeros) |
| 14 | `34`@0 -> `34`@2 | 2 moves (rule `[2,2]`, then rule `[-7,4,-7]`): the a = 4 case is shorter, consistent with E-067's `34 -> 44` |
| 13 | `5046`@0 -> `5046`@2 | 7 moves: 50460000000, 40360003000, 40350033300, 40340334000, 40333350000, 33340350000, 03500350000, 00504600000 |
| 12 | `4`@0 -> `3`@1 | one rule move (0300000000) |
| 12 | `5`@0 -> `3`@2 | orbit of `5`@0 closed at 148 rows, target absent |

Single-step neighbourhoods (`theorist_k0moves.py 12 3 7`): `3@0` 2 neighbours (orbit 10), `4@0` 5 (1410), `5@0` 2 (148), `6@0` 2 (492), `7@0` 2 (158). Orbit statistics of P, Q at n = 12..14 (`theorist_orbitstats.py`) show no separating invariant.
Weakest points: the path is a shortest path of the move set used (rules + edge + double + spectator), not a hand proof of each rewrite; "a-1 steps" is read from five values of a, n = 12..17.

## Reproduction

```
.venv/bin/python workshop/rounds/013/theorist_path.py 12 35 0 35 2      # also 12 35 1 35 3; 13 36 0 36 2; 15 37 0 37 2; 16 38 0 38 2; 17 39 0 39 2; 13 5046 0 5046 2  (each < 10 s, n = 13 5046: 1 min)
.venv/bin/python workshop/rounds/013/theorist_path.py 12 4 0 3 1        # 1 step
.venv/bin/python workshop/rounds/013/theorist_path.py 12 5 0 3 2        # closed, absent
.venv/bin/python workshop/rounds/013/theorist_k0moves.py 12 3 7         # neighbours of k@0
.venv/bin/python workshop/rounds/013/theorist_orbitstats.py 12          # also 13, 14
timeout 10m .venv/bin/python workshop/rounds/011/theorist_word.py 17 5046   # about 3.5 min; same for 5056
.venv/bin/python workshop/rounds/013/theorist_size.py 15 5046 0 100000      # sizing (19 s)
```

## Prior record

E-079 lists the orbits P, Q and says no move sequence was exhibited; E-067 has `34 -> 44` and drift families; E-070 says the parity translation belongs to `5046`/`5056` at odd n (grep "5046" in EXPERIMENTS.md).
Not found in `research/`: the staircase paths, the "only width-4 rules have a lone-relation left side" account of `4@0 -> 3@1`, or `5046`/`5056` at n = 17 (E-076 had "nothing for n >= 17"). Nothing in `RETRACTIONS.md` involved (grep of "5046": none).

## Code changed

New scripts only: `workshop/rounds/013/theorist_path.py` (labelled BFS), `theorist_k0moves.py`, `theorist_orbitstats.py`, `theorist_size.py`. No library change, no tests run.

## Next

- Skeptic: re-run two paths and check, from `lnaMoves` matching, that each double move `[t,t]` is applicable as stated; test the template for a = 10, 11 (n = 18, 19) and the claim "a-1 steps".
- Theorist: try to prove that offset shift by 1 is impossible via an invariant of the double-mutation moves; the paths suggest looking at the multiset of mutated vertices (each double move mutates a vertex twice).
- Experimentalist: `orbitclass` for n = 17 can be replaced by running only the 5046/5056 orbits (done) plus the census-rule words; the remaining 137 cores are the open part. No overnight run needed for A/B.
