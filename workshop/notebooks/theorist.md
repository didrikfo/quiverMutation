# Theorist notebook (rewritten round 013)

## What I believe now
- P (orbit of `3@2`/`3@3`) and Q (`5@0`/`6@0`) are each closed under offset shift by 2. Explicit: `3a`@o -> `3a`@(o+2) in a-1 double-mutation moves, a = 5..9 (`theorist_path.py`, rounds/013);
  `35`@0 -> `35`@2 at n = 12 is 4 moves (+2 spectator; rule [2,2]; doubles [5,5],[4,4],[3,3]). Shift by 1 unreachable (orbit closed) but no invariant or proof found: "parity class" names two orbits, not an invariant.
- `4@0 -> 3@1` is ONE move: an anchored width-4 rule in `lnaMoves.ALL_MOVES` (4 -> 3). Only width-4 rules have a lone-relation LHS; so `k@0` (k >= 5) has only 2 neighbours (`kk`, `k0..03`). Explains letter 4, not `444`.
- `5046`, `5056` at n = 13, 15, 17 also have two parity orbits (n = 17: 122673 / 54266, one key; sizes differ so mirror cannot merge) and a staircase path at n = 13; they sit in no single-relation orbit. Growth about 4.3 per n -> n+2 (orbit walk 3.5 min at 17).
- Earlier (round 011): census rule reproduces A and B minus `5046 5056`; `406` excluded by mirror clause; nulls: no GF(2) functional, integer statistic or SNF separates P from Q (do not retry). Rule was fitted at n = 12..16 (E-077).
- Round 009: drift family `aax`; a = 4 special (`444 -> 34 -> 44`); `34`@0 -> `34`@2 is only 2 moves (shorter than the a >= 5 staircase).

## What I tried
- Round 013 scripts: theorist_{path,k0moves,orbitstats,size}.py. Labelled BFS = quickest way to see mechanisms; orbit stats (nrel, max/min letter, sums) show nothing.
- Not done: proof of "no shift by 1"; a = 10+ ; other 137 cores at n = 17; mirror membership of the n = 17 5046 orbits stated only by size argument.

## Next
- Look at the multiset of mutated vertices / mutation-sequence parity along the staircase (each step mutates one vertex twice) for an invariant; try a signed count of mutations per vertex mod 2 over all moves.
- Compare P/Q path template to the 7 of E-059 (`344 366 4044 ...`): do they get a staircase too?
- Blind spots: "a-1 steps" is read off 5 values; shortest path in my move set is not a proof of derived equivalence step by step (moves are table-verified, not re-derived by me).
