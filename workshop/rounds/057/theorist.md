# T1: s = n - k(c) holds on the cores checked, no linear letter rule for k(c) found; T2: the 333@0 row set is the 444 orbit, and it is huge

author: theorist · round: 057 · kind: proposal
thread: T1, T2 · bears on: H-021, E-061..E-064, E-067, E-090
scope: breadth slot, no new census. Linear fits on 17 recorded k values (11 independent of the 33x family); row-set check n = 13..16 with full tables; word breakdown of the 444 orbit at n = 14 only. Orbits of the reduced walk only, not derived equivalence.

## Response to referee

1. Done. Prior record now cites E-090 (J labels, `theorist_label.py`), E-067 (`444@o` in the orbit of `333@0`) and E-088/E-093; "not found" dropped. The row-set equality is E-090's J-label made explicit; only the set comparison is new.
2. Done. Full tables for n = 13, 14, 15, 16 below (`theorist_same_n{13,15,16}.txt`). The other `333` placements are disjoint from the `444` orbit at all four n.
3. Done. Refit without the 33x block (11 cores, `theorist_kfit.py no33`): best leave-one-out SSE 48.3 (`last`, `lastm1`), at most 3 of 11 exact; so worse than on 17, as expected once the linear 33x block is removed. Ledger wording changed to "no linear letter-statistic fit predicts k(c) on 11 independent values; a drift-based rule untried", not "no rule".
4. Source: round 004 `experimentalist_shift_table.txt` has n = 16 rows (nofit for 3355, 3445) and E-064 gives `46` nofit at 16 (9 of 12). Cited here, not rerun.
5. Done at n = 14 (`theorist_breakdown.py 14`, `theorist_breakdown_n14.txt`): the 3767 rows have 2871 distinct trimmed words and 3442 rows have a trimmed word of length >= 5. The `33y` rows are about 20 (333 twice, 334..339 once each, ...). So the orbit is not "33y plus a few neighbours"; the upper bound of E-067 concerns only the `33y` rows and the orbit is far larger. The proposed "remainder is 44y, 345, 3334-type" is refuted. Other n not run.

## Claim

**T1 (H-021', s = n - k(c)): propose closing it for the next ledger** (chair's call), as "n-independent on the 12 cores checked (n = 13..15, 9 of 12 at 16); no linear letter-statistic fit predicts k(c)". The only thing that would have kept it open is a rule for k(c) that predicts a held-out core (maverick, round 047). I tried the cheapest family (k linear in 1-3 of: length, digit sum, first, last, max, min nonzero, zero count, last gap) with leave-one-out on 17 cores: best squared error 54, at most 6 of 17 exact; on the 11 cores outside the 33x family best 48.3, at most 3 of 11 exact. No cheap rule in these features. Not claimed: that no rule exists.

**T2 (centre s(c), k(33x) = 2x): do not close; one concrete reopening question, now partly answered.** The recorded flag was "33x orbit size 3767 vs the n = 14 `444` orbit". Settled: the orbit of `333@0` and of `333@(n-6)` is the same row set as the `444` orbit (the orbit of `444` at its middle offset), at n = 13, 14, 15, 16 (sizes 2386, 3767, 5648, 8134). The other `333` placements are disjoint from it. So the c = 3 chain of E-067 is the `444` orbit, which explains why `44x` "attaches to the c = 3 orbit" (E-067 (4)) and why the argument for `44x` cannot give k = 2x - 1. Consequence for the ledger: E-067's upper bound is false as a statement about the orbit, true only for the `33y` rows in it. The sentence "k(33x) = 2x" is a statement about which `33y` rows share an orbit, not an orbit description.

**Breakdown done (n = 14):** the orbit has 2871 distinct trimmed words in 3767 rows, 3442 rows of length >= 5; the `33y` rows are about 20. So E-067's claim is only about which `33y` rows share an orbit, and can be closed as such ("lower bound derived, `33y` rows exhaustive at n = 14..17"). The remaining open question is the one the orbit size poses, not T2's: what structure makes the `444` orbit so large (H-020 / T3 territory). Proposal: close T2 as stated.

## Evidence

T1 fit (`theorist_kfit.py`, seconds): data are E-064's shift table at n = 13 (12 cores: 334 3344 3345 3355 3444 3445 4045 45 4506 4556 46 504) plus `33x` x = 3..8 from E-067 (k = 2x; 334 and 333 appear in both; the 33x rows are counted once). Models `k = a + sum b_i f_i` with 1, 2 or 3 features; best two-feature fit `first, last` has leave-one-out SSE 54.4 and predicts 5 of 17 exactly. Single feature `last`: SSE 58. Compare: the spread of k is 6..16, so these fits barely beat the mean. The mechanism that does give k, the drift to a self-dual seed (E-067), is not a function of simple letter statistics: `44x` shares the drift and fails.

T2 check (`theorist_same.py N`): per placement of `333`, orbit under `freeMoves.REDUCED` (limit 300 000, all closed), compared as a set of rows with the `444` orbit.

| n | |444 orbit| | 333@o equal to it | other 333 placements |
|---|---|---|---|
| 13 | 2386 | o = 0, 7 | o = 1..6 sizes 447 449 501 501 449 447, overlap 0 |
| 14 | 3767 | o = 0, 8 | o = 1..7 sizes 320 886 364 491 364 886 320, overlap 0 |
| 15 | 5648 | o = 0, 9 | o = 1..8 sizes 763 743 881 877 877 881 743 763, overlap 0 |
| 16 | 8134 | o = 0, 10 | o = 1..9 sizes 516 1416 594 1640 310 1640 594 1416 516, overlap 0 |

(Sizes of the 444 orbit match the E-093 `|S|` series.)

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/057/theorist_kfit.py [no33]   # seconds
timeout 10m .venv/bin/python workshop/rounds/057/theorist_breakdown.py 14   # about 1 min
timeout 10m .venv/bin/python workshop/rounds/057/theorist_same.py 14  # about 40 s; 13, 15, 16 took about 2 min together
```

## Prior record

E-090 labels each closed orbit by the `333` offsets it holds (`J = {0, n-6}` for `444`, n = 12..17); E-067 (line ~858) says `444@o` lies in the orbit of `333@0`; E-088/E-093 give the sizes. So the equality is recorded; new is only the explicit row-set comparison at n = 13..16 and the n = 14 word breakdown. The T1 closure follows maverick's round-047 note and the round-048 header rewrite of H-021; the failed fit is new but weak (17 points, linear models).

## Code changed

None. New scripts: `workshop/rounds/057/theorist_kfit.py`, `theorist_same.py`.

## Next

- Chair: ledger line for T1 "closed: no rule for k(c)" if the skeptic agrees the fit is a fair attempt; T2 stays open on the question above.
- Weakest point of the T1 closure: only linear letter statistics were tried; a rule from the drift structure (which seed, how many steps) was not. If anyone wants to reopen T1, that is the way, and it is the E-067 programme for `34x`/`345`.
- T2 breakdown at n = 15, 16 not run (n = 14 done).
