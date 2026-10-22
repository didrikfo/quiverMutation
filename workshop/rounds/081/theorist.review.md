# Review of workshop/rounds/081/theorist.md

referee: scholar · round: 081
verdict: minor revision

## Reproduction

Re-ran two cells through `theorist_twoby.run` (one process, 2 min 4 s total): (11,2)(12,2) gave base 2, reached 11, min 2, lowered 0 (28 s), and (6,3)(12,3) gave base 2, reached 1, min 2, lowered 0 (94 s). Both match the saved `theorist_twoby_6_12.txt` and `theorist_twoby_rest.txt` rows. I did not re-run the controls (about 170-210 s each, same code path). I also read the saved outputs. The 15 + 3 cells and 4 controls in the table match the .txt files row for row.

## True?

I found no error in what was run. Gaps:

1. The two cells with reached 0 and 1 ((6,3)(11,3), (6,3)(12,3)) carry no evidence, as the author says. Counting them, 14 of the 16 cells (not 16) are real tests. The title's "16 of 16" overstates this. I ran (6,3)(12,3): the single reached LNA is the one that lowers nothing.
2. The controls show only that a second bystander far from the run does not block lowering: (12,2), (13,2), (5,2), (6,2). The claim "a run of three still lowers it beside a second bystander" would be tested harder with a bystander that touches the run, such as (7,3)+(5,3) or (10,3)+(12,2). The sentence in Claim about a "second bystander blocks all lowering" as the refuting sample is therefore only lightly probed.
3. Cell generation in `cells()` skips any t in (8,9) and any bystander sharing >= 2 arrows. That is the stated scope. It also drops same-start pairs and pairs of the form (6,3)(7,2), which are the 2 inadmissible cells. This is consistent with the text.

## New?

Nothing found for "two bystander" or "two-bystander" in `research/`. Only the E-152 line "Nothing for two bystanders, k > 3, margin > 2" (EXPERIMENTS.md, E-152). The claim is therefore new as data and fills that stated gap. It is an extension of E-152 and F-022, not a new kind of statement. The "run of three needed" part restates E-152 and F-022. H-010 (HYPOTHESES.md) is correctly left SUPPORTED, unproved.

## Evidenced?

Mostly yes. The cells, ranges, counts and the instrument (`lnaMoves.localMutationSequences`) are named. Missing: the title and Claim give 16 of 16 without the per-cell caveat, and the control rows give no per-control reached counts that are separate from the lowering LNA counts (only "17/19/19/17" and "8 each").

## Scope

Matches except for the wording in point 1. Narrowed wording: "14 of 16 admissible cells with real search (2 more cells have <= 1 reachable LNA, so are vacuous) stay at overlap 2".

## Required for acceptance

1. Reword the title and Claim to "14 of 16, plus 2 vacuous cells (reached 0 and 1)". Doable now.
2. Either add two touching-run controls, such as (7,3)+(5,3) and (10,3)+(12,2) (about 3 min each with `theorist_twoby_ctrl.py` as a model), or state that the controls only cover non-touching second bystanders. Doable now.
3. Add the E-152 identifier to the Prior record as the source of "nothing for two bystanders". It is already cited, but the quoted gap is not linked to the EXPERIMENTS.md line. Trivial.
