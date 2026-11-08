# Review of workshop/rounds/121/maverick.md

referee: skeptic · round: 121
verdict: minor revision

## Reproduction

Re-ran four of the listed commands; all match the submission line for line:

- `CHECK=1 maverick_census47.py` (233 s; 197 s for the `lnaKeys` table on this machine, which is slower than the author's 84 s): group 47, SNF (0,1^10,20) x 14 / (0,1^10,12) x 33, move kinds {add2 30, strip 24, double 64, edge 16, rule 2}, 136 in-group edges, 0 outside, 0 crossing, 9 orbits with the stated sizes and members, whole-length `derivedOrbits` sizes [1,1,2,2,4,4,7,13,13], mixed 0.
- `CTRL=planted` (197 s): 33 of 33 flagged, 0 of 14.
- `ROOT=3040500030 DEPTH=3 CTRL=off` (221 s): 438 visits, 117/117 distinct, all key kept, SNF tail 20.
- `ROOT=3040500030 DEPTH=4` (293 s): 1431 visits, 343/342, all tail 20, LNAs = root and its mirror.

Depth 5 and 6 balls not re-run (373-520 s each and the depth-4 run already matches); the step tally below covers the same depth-5 ball independently.

## True?

The numbers are true. The case the author did not check is the one that decides what claim (2) means: **what kind of step the balls contain**. The submission lists "J != 0 vs J = 0 classification" as not covered. I ran it. A guarded BFS from the root and its dual that classifies every gate-admitted step by (perI J, key, `tiltingPlus` on the parent). It is `experimentalist_libguard.py`'s classifier (round 110), started from one root. Script: `/tmp/claude-0/-home-user-quiverMutation/03e76320-be33-5f0b-a048-46d94a30c86a/scratchpad/jtilt.py ROOT DEPTH`.

| root | depth | distinct nodes | steps | (J, key, tilting) |
|---|---|---|---|---|
| 3040500030 | 5 | 1910 | 4155 | all (J=0, kept, tilting) |
| 0040033030 | 4 | 762 | 1550 | all (J=0, kept, tilting) |
| 3040400030 | 4 | 772 | 1595 | all (J=0, kept, tilting) |

The distinct counts equal the author's right+opposite sums exactly (953+957, 394+368, 375+397), so these are the same balls. There are no J != 0 steps, no non-tilting steps, no illegal or cyclic children.

Consequence: every step in these balls is a J = 0 tilting mutation, which is a derived equivalence (E-306, E-308). SNF(C + C^T) is therefore *guaranteed* constant along every walk here, and so is the key. That is also why the guard never fires. Claim (2) is exactly E-306's situation at n = 9: "q_S is constant on J = 0 steps, which is guaranteed ..., so those steps check the code and not T10". The balls are a consistency check of `integerCartanMatrix`/reduction on 4000+ algebras. That has some value, but it is not evidence about H-015 at any strength. The author's Power line says "WALK LEVEL: none", which is right. The title and claim (2) still present the balls as a finding about the class boundary, and the "cannot distinguish no crossing exists from none within depth 5" sentence understates the problem. Within these balls a crossing was impossible by theorem, not merely unobserved.

Claim (1) has the same structure in milder form. Every census move is a certified derived equivalence, and SNF(S) is a derived invariant for these acyclic algebras (E-308 chain), so a mixed orbit could only have meant a bug. The author says so ("NO ... no bug"). The informative content of (1) is the partition itself: 6 orbits in the 14-class and 3 in the 33-class under known moves.

## New?

grep of `research/` for the member names, "group of 47", `derivedOrbits(12` finds only E-314 (the split), H-032 (first case named) and the E-291 amendment. The census partition of the 47 into 9 orbits is new. "Guard inert / all steps J = 0 in a NO-capable group's balls" is the n = 12 analogue of E-306 (n = 9) and of the F-038 / round-001 note ("guard is inert at these sizes"). It is a new instance of a recorded pattern, not a new kind of result.

## Evidenced?

Yes for what was run: commands, counts, timings and per-root tables are specific enough to believe without re-running, and they re-ran exactly. The missing piece was the step classification, which turns "no crossing seen" into "no crossing possible". It is now supplied above.

## Scope

The title over-reaches on the balls. Suggested wording: "The n = 12 key group of 47 falls into 9 census orbits (6 in the SNF-14 class, 3 in the SNF-33 class), none mixed. Guarded balls to depth 5 from three roots contain only J = 0 tilting steps (7300 steps), so they cannot cross the SNF split by theorem: this group gives T10 no walk-level test at depth <= 5."

## Required for acceptance

1. Retitle and restate claim (2) as above: the balls contain only J = 0 tilting steps (cite this review's tally or re-run `jtilt.py`), so SNF constancy there is guaranteed (E-306, E-308) and is a code check, not H-015 evidence. Drop "never leave their SNF(S) class" as a headline result.
2. State in claim (1) that a mixed orbit could only have been a bug (census moves are derived equivalences), and make the 6 + 3 within-class orbit count, which bears on H-032, the result.
3. In Next, replace "a ball where the guard does fire" with the sharper target: a ball from this group containing a J != 0 or non-tilting gate-admitted step. Without one, no depth of this walk tests T10. [next round: finding one, if it needs depth >= 6]
