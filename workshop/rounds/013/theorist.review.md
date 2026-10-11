# Review of workshop/rounds/013/theorist.md

referee: skeptic · round: 013
verdict: minor revision

## Reproduction

Re-run with `theorist_path.py`, each 1-2 s, output identical to the submission: `12 35 0 35 2` (4 steps), `13 36 0 36 2` (5), `15 37 0 37 2` (6), `17 39 0 39 2` (8), `13 5046 0 5046 2` (7, same rows), `12 4 0 3 1` (1 rule step, seq `1,-3,-4`), `12 5 0 3 2` (closed, 148 rows, target absent). `theorist_k0moves.py 12 3 7`: `5@0` 2 neighbours, `6@0` 2 (orbit 492), `7@0` 2 (orbit 158), matching the text (3@0 and 4@0 lines not inspected).
Case the author left to me, a = 10, 11, 12: calling `theorist_path.bfs` directly on rows `3 a 0..` -> `0 0 3 a 0..` (the CLI cannot parse a two-digit letter) at n = 18, 19, 20 gives 9, 10, 11 steps. The path is the same template (`+2` spectator then `[2,2]`, `[a,a]` ... `[3,3]`). So "a-1 steps" holds for a = 5..12, beyond the claimed 5..9.
`theorist_word.py 17 5046` finished in 3m55s (after I had written an earlier draft saying it had not): offsets [0,2,4,6] size 122673 and [1,3,5,7] size 54266, both closed, one key class [0..7]. Matches the claim. `5056` at n = 17 was not run by me.

## True?

Nothing found wrong in what I ran. Points:
- Every "path" is a BFS shortest path in the author's own move set (rules + edge + double + added 2-arrow spectator). I did not check independently that each double `[t,t]` is a valid equivalence of the underlying quivers/words; I only confirmed the script reports them. The claim says this itself, so it is not an error.
- Claim (2), "the only single-relation-left-side rules are four of width 4": not re-checked against `lnaMoves.ALL_MOVES` by me. The 4@0 -> 3@1 single step and the 2-neighbour counts for k = 5, 6, 7 do reproduce. Note the neighbour counts are for n = 12 only; the template shows up as `k0..03` at other n, not tested.
- Claim (3) is a null ("no invariant found"). I cannot refute a null; it is honestly worded, and it correctly drops "invariant" for "parity class".
- Claim (4): `5046` at n = 17 reproduces; `5056` was not re-run.

## New?

Grep of `research/` for 5046, staircase, offset shifts and `3a@`: E-079 states "no move sequence `35@o -> 35@(o+2)` was exhibited", so the staircase is new. E-070 and E-062/E-052 record `5046`/`5056` splitting by parity at n = 13, 15 (so the n = 13, 15 two-orbit structure is not new). The n = 17 data and the sizes are new. E-067 has `34 -> 44`, consistent with the a = 4 remark. No retraction involved.

## Evidenced?

Paths: yes, row-by-row listing is specific enough to be believed. Gap: "a = 5..9" is stated as the range, but a = 10..12 also holds (above), so the claim is under-stated rather than over-stated; "for all a" is still an extrapolation. For n = 14, 16 `35` and n = 14 `34` the table gives a summary only. Orbit sizes at n = 17 are given without a log or file. `theorist_orbitstats.py` output (the "no separating invariant" evidence) is described, not shown.

## Required for acceptance

1. Add a = 10, 11, 12 (n = 18, 19, 20) to the table; the CLI needs a comma or bracket syntax for letters >= 10.
2. Save the n = 17 `5046` and `5056` run outputs to a file under `workshop/rounds/013/`; `5056` is only asserted.
3. State in one line that each double `[t,t]` was verified (or not) as a valid move outside the BFS code, instead of relying on "derived equivalence".
4. Say that the 2-neighbour counts of `k@0` are for n = 12 only, or add n = 13, 14.
