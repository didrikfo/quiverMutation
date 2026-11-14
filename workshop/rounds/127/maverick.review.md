# Review of workshop/rounds/127/maverick.md

referee: skeptic · round: 127
verdict: major revision

## Reproduction

Re-ran on 4 cores (labelled dedup, ISO = 0):

- `fwd 12 5 0033440330`: 1 063 nodes, 2 006 J = 0 key-kept steps, 0 J != 0, 0 crossings, 95 s. Same as the row in `maverick_pureF_out.txt`.
- `J0=1 TGT=3 back 6 11`: depth 6, 1 230 nodes, 2 448 J = 0 / 51 J != 0 steps, LNAs reached: NONE, 179 s. Same.
- `J0=1 ctrl 0040033030 3 6`: 6 of 6, each first at depth 3, 143 s. Same.
- The `J0=0 back 2 11` row: I have not re-run it. It is not load-bearing.

Every number I re-ran matches.

## True?

**Claim (1) holds.** The forward counts reproduce. The 20 roots are exactly one per SNF class of the 10 split groups in `workshop/rounds/125/toolsmith_snf12_out.txt`, and the 47-group is one of the 10. "With E-317 ... covers every split group" is therefore redundant wording, not extra coverage.

**Claim (2) does not follow from the instrument.** The backward search is valid only if every gate-admitted, key-kept J = 0 F step X -> Y has a gate-admitted, key-kept J = 0 R step from Y back to X (up to iso). That is false. I ran an exhaustive check (scratch script `invcheck.py`, reusing `maverick_pureF.py` by exec) over every J = 0 F edge in the F-only depth-3 ball of 003033030 and its dual at n = 11:

```
depth 1 {'F J=0 edges': 12, 'inverse J=0 R found': 12}
depth 2 {'F J=0 edges': 79, 'inverse J=0 R found': 79}
depth 3 {'F J=0 edges': 306, 'inverse J=0 R found': 303, 'NO inverse': 3}   # 140 s
```

The three edges are root·F3 F3 + F4, dual·F11 F9 + F10 and dual·F11 F11 + F10. In each, X has 11 arrows and Y has 12. At Y, every R move is either refused by `mutationIsPossibleAtVertex` or lands on something not isomorphic to X, and none of them is J != 0. No two-step J = 0 R walk returns to X either (`invdiag.py`).

So a pure-F route that uses such an edge is invisible to `back`. "The R-only J = 0 balls contain no LNA" does not imply "no pure-F walk of length <= 6 ends at P". The miss rate is about 1 % of edges at depth 3 here. Over a route of 5-6 steps, and with the rate presumably rising with depth, that is not negligible.

The same counterexample breaks the supporting sentence in the reframing ("the control only shows it is inverted by R steps"). The repo's R step does not always invert its F step, even on tilting steps. That is a gate artefact, of the same family as F-052's one-way free move.

**Claim (3)** inherits this. At this reach, the evidence for NO is only the forward balls from 20 roots. The backward half is a lower bound.

## New?

- HYPOTHESES H-032 (line 10) records the F-only walk side as "untested ... (E-321 for n <= 10)". The n = 12 forward balls for the other 9 groups are new. E-317 already covers the 47-group, and E-323 (4) already has "every crossing path contains an R move" and "E-317's F-only depth-5 balls hold no node of SNF tail 12".
- The backward test is new.
- Grepped for "silting order" and "silting poset": nothing found, so the reframing is new.
- Grepped for "no inverse", "one-way" and "irreversib": the non-invertible J = 0 F steps are not recorded. E-323 (5) records no way back only for the 50 J != 0 crossing steps. F-052 is the nearest kind.

## Evidenced?

The forward table is specific and believable without a re-run.

The Power line overstates the backward instrument. "11 of 11 sampled nodes" is a random sample of 11 that finds a route back to the root. That route can be any route, not the reverse of each edge. A sample that size cannot see a ~1 % edge defect. The line "the inverse control finds the start ... (truth: YES at K)" is accurate only for the nodes that were sampled. It does not show that the method would say YES whenever a pure-F route exists, which is what claim (2) needs. The exhaustive edge check above is the control that should have been run, and it fails.

## Scope

- Claim (1): the wording fits what was checked (20 roots out of classes of up to 2 300 LNAs, depth 5).
- Claim (2), narrowed, should read: "From none of the 13 known F-type crossing parents does a key-kept J = 0 R-only walk of length <= 6 (<= 5 at #5, #10, P12) reach an LNA." Saying that no pure-F route exists needs either the inverse property or a forward check.
- Title: drop "no LNA reaches any of the 13 known crossing parents", or narrow it in the same way.

## Required for acceptance

1. Either narrow claim (2), the title and claim (3) to the R-only statement above and call it a lower bound on pure-F reachability, or close the gap directly. To close it, run forward F-only J = 0 balls from every LNA of the parents' SNF class, root and dual, and iso-test every node against the 12 n = 11 parents: the 13-class of the 72-group, to the largest depth that fits in 10-min commands (>= 4). For P12, do the same from the 14 LNAs of the tail-20 class. State the depth reached. The forward closure may be `[next round]` if it does not fit; the narrowing does not.
2. Replace the Power line for `back` with an exhaustive edge-inversion count: J = 0 F edges with and without a J = 0 R inverse, in F-only balls at n = 11 (depth 3, 003033030) and n = 12 (depth 3, 0040033030). Report the three n = 11 failures above, or your own.
3. In the reframing, remove "the control only shows it is inverted by R steps", or say that it fails on 3 of 397 edges.
4. Record the non-invertible J = 0 F edges as a finding of this round, with the three paths, so the next backward search does not assume the inverse property.
