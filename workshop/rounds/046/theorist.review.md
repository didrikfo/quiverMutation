# Review of workshop/rounds/046/theorist.md

referee: Maverick · round: 046
verdict: accept, narrowed (observations only; H-010 stays unproved)

## Reproduction

- `theorist_single.py 3 2 10:4`: 24 s. `[(8,3),(9,3),(10,4)]`, steps=3, margin=2, base=2, reached=4, min=0, lowered=1. Matches the author. So `(10:4)` lowers at k = 3 and the k = 2 non-lowering is not a counterexample to "a run of three lowers by k = 3".
- `theorist_locality.py 4 8` (n = 4..8, under 10 s): events 20, 70, 252, 924, 3432; LNA results 10, 28, 84, 264, 858; offsets only -2, -1, 0, ratio 1:2:1 (n = 8: 297 / 594 / 297). Matches `theorist_locality_4_10.txt`. n = 9, 10 not re-run (n = 10 is 3.5 min by the author's count); the pattern is identical in shape, so no smaller case was needed to expose an error.
- The k <= 2 table was reproduced by me in round 045 and is unchanged. The left sweep and the other k = 3 cells were not re-run; I take them from the committed outputs (the numbers agree with the 045 variation I ran: `(7:3)` 4 reached / 2 lowered / min 1 at k = 2).

## True?

Yes, as far as it was run. Points checked against my 045 list:
1. 21 of 21 is correct (23 rows minus the 2 runs of three).
2. Per-cell placement counts at k = 3 are stated, and the cap-cut cells are marked "not run".
3. The `(10:4)` result is reported as lowering, and the claim is correctly reworded to "necessary in every cell tested, sufficient in 5 of 5 tried at k = 3". Not "criterion" any more.
4. Margin 2 is in the Scope and the Claim.
5. Locality is marked a conjecture for k steps; L1 (one step, LNA to LNA) was run and is correctly said not to cover non-LNA intermediates.

Remaining small points, none blocking:
- The L1 "changed start" test is about relation starts only; I did not check that it also covers a changed relation's length or the appearance/disappearance of a relation. The author says "changed relation starts"; the claim is worded that way, so it stands as stated.
- "5 of 5 runs of three tried" includes `(8:4)(9:4)(11:3)` where the pair itself is `(8:4)(9:4)`, a different shape. It is listed separately; fine.
- Item 7 withdraws the E-068 sentence honestly. The F-025 "cheaper line" is still not engaged, said so.

## New?

- Bystander result: F-022's table already holds the statement (`(1:3)(2:3)(3:3)` to 0; `(4:2)`, `(5:2)` inert; `(1:4)(2:4)(3:4)` only to 2); E-025 qualifies it for long runs. The author now says so (item 3). What is new is modest: the sweep over gap, length, side and three pair shapes, and `(10:4)` needing 3 mutations. Not a lemma.
- L1: grep of `research/` for locality, "changed start", "one-step", "exactly two" finds no count of the kind. Nearest priors are E-021/F-022-era `probe.py` locality remarks (FINDINGS around F-022 end: "locality is what a move rule is"), and the square-structure table (EXPERIMENTS ~3180, FINDINGS ~1716: a mutation leaves the line only into a square with short side exactly two; 336 "a line again"). Those are about rules and squares, not a count over all LNAs of one-step LNA-to-LNA results. The count (exactly 2 of n mutations per LNA stay LNAs; offsets 1:2:1) is a new, checkable, small observation. It is probably explicable by the square picture (a line returns to a line only at a vertex adjacent to the relation end), which the author did not link; I suggest the link, not require it.

## Evidenced?

Yes for what is claimed. Outputs and scripts are committed, commands are listed with times, and the claims are bounded by the table. The only unevidenced piece is the cap-cut k = 3 rows, honestly marked.

## Scope

Title and Claim now match: k <= 2 full, k = 3 partial, margin 2, one bystander, pair `(8:3)(9:3)` plus sampled shapes, L1 for n <= 10. No narrowing beyond one wording: "the same holds for the other pair shapes sampled" should read "consistent with" because the other shapes were run at k = 3 only, with one lowering example (`(11:3)` on `(8:4)(9:4)`).

## Required for acceptance

None. Promotable (to FINDINGS/EXPERIMENTS as an observation, level "tested on small cases"):
1. For the pair `(8:3)(9:3)` with >= 6 empty vertices each side, margin 2, one bystander m = 2..4: no bystander sharing <= 1 arrow lowers the overlap (21 right and 14 left placements at k <= 2; 10 right placements and 2 left at k = 3); runs of three lowered by k = 3 in 5 of 5 tried, `(10:4)` only at k = 3 (not k = 2). A refinement of F-022, not new in kind.
2. L1: for all LNAs n = 4..10, exactly 2 of the n single mutations give an LNA, and the changed relation starts lie in {v-2, v-1, v} (counts per n above). Computation, not a proof; says nothing about k-step windows.
3. Status: H-010 unproved, no counterexample. Not promotable: any claim of a k-step locality window, anything for two bystanders, k > 3, margin > 2, or long runs (E-025).
[next round] cap-cut k = 3 rows; two bystanders; explain the 2-per-LNA count via the square picture.
