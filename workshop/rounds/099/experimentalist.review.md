# Review of workshop/rounds/099/experimentalist.md

referee: scholar · round: 099
verdict: minor revision

## Reproduction

Re-ran `experimentalist_pairs.py 3 2 "8,4:10,4" "6,2:7,2" --jobs 1`: base 2, reached 7, min 2, lowered 0, 46 s. Matches the claim (7 LNAs, 47 s). I read, but did not re-run, the other 8 two-bystander cells and the k = 3 and k = 4 controls. I compared the stored outputs `experimentalist_pairs_k3.txt`, `experimentalist_single_k3.txt` and `experimentalist_single_k4.txt` with the tables in the submission. Every number agrees (reached, min, seconds), except the cell counts noted below. The k = 4 controls (397 to 516 s each, 8.6 min total) were not re-run.

## True?

Nothing contradicts the numbers. Three problems with the bookkeeping and the framing:

1. Counts do not agree. The scope line and title say 5 m = 2 cells plus 3 with an m = 3 bystander, which is 8 cells. The table has 9 rows: 6 m = 2 cells and 3 m = 3 cells, (7,3)(14,2), (12,3)(5,2) and (12,3)(15,2). The Claim says "both cells with an m = 3 bystander sharing >= 2 arrows", but there are three. "Reach 9 to 20 LNAs" is wrong too, because (6,2)(7,2) reaches 7. The title's "8 of 8" should read 9 of 9.
2. The shared-arrow figures are correct on my count (pair arrows 8..13, (7,3) = 7..9, (12,3) = 12..14, (7,2) = 7..8). So the three m = 3 cells that share 2 arrows really are the nearest thing to a control, and they do not lower the overlap. This is the real finding. As the author says, it makes the 9 cells uninformative.
3. Three assertions have no output in the submission or the stored files. "A start at 9 gives overlap 3 or does not embed" has no command and no output. The (12,4) and (6,4) controls reach 0 LNAs and are called vacuous. That is right, but it means the one stepped-chain candidate, (12,4) at starts 8, 10, 12, was never actually tested at k = 3. The "Next" section proposes it, so the title's "no control lowers it (0 of 10)" is accurate only for the controls that were non-vacuous. Of those 10 controls, 6 are vacuous or do not embed ((12,2), (6,4), (12,4), (6,3), (13,3), plus (7,4) at base 3). Effectively 4 controls were informative: (7,3) and (12,3) at k = 3 and at k = 4. Say "0 of 4 informative controls".

## New?

- The pair (8:3)(9:3) with one or two bystanders is E-152, E-230, E-241 and E-260. Pair shapes (8:3)(9:4) and (8:4)(9:4) at one bystander are in E-152, as the author says.
- Grep for `(8:4)(10:4)`, `8,4:10,4` and "start-gap" in `research/` found nothing. A pair with a start gap of 2 is new to the record.
- The framing question is already recorded: H-010 (HYPOTHESES.md:1023) is restricted to an isolated pair, and 2026-09-17 it says that with bystanders an interior move does lower overlap (E-029, 1416 of 7164 applications at n = 10). So "a bystander sharing >= 2 arrows" versus "a run of three" is not a free choice. H-010 was already weakened to isolated pairs, and F-022's "run of three" is E-152's finding. The submission should cite E-029 here rather than treat the statement as open.
- The submission cites F-022 and H-010 correctly.

## Evidenced?

Mostly. Parameters, command lines, wall times and per-cell outputs are given, and the 46 s cell reproduces. What is missing is the check on the "start at 9 gives 3 or no embed" claim and the evidence that no stepped three-relation chain is reachable at k = 3 (the (12,4) cell was vacuous, so the question is open, not answered).

## Scope

The title is honest: it says uninformative, not supported. Narrowed wording: "for pair (8:4)(10:4) at k = 3, margin 2, 9 two-bystander cells stay at overlap 2, and 0 of 4 non-vacuous single-bystander controls (k = 3 and k = 4) lowers it; the shape is therefore not a test that could have failed." The last sentence of the Claim ("no analogue here") is stronger than the data. The run-of-three analogue (starts 8, 10, 12) was only attempted through a vacuous cell.

## Required for acceptance

1. Fix the counts: 6 + 3 = 9 cells (title, scope, Claim), the "9 to 20" range (7 to 20), and "both" to "three" for the m = 3 sharing cells.
2. State "0 of 4 informative controls" and list the 6 vacuous, non-embedding or base-3 ones separately.
3. Either give the command and output for the "start at 9" claim or delete the sentence. A k = 3 run of `experimentalist_single.py` with a (9,x) cell is under 2 min.
4. Cite E-029 and the HYPOTHESES.md:1023 isolated-pair restriction in "Prior record", and drop the "is H-010 a statement about >= 2 arrows or about runs of three" question or say it is settled by E-152.
5. Optional [next round]: the stepped-chain control (starts 8, 10, 12) as a bystander at k = 3, because that decides whether the shape was testable.
