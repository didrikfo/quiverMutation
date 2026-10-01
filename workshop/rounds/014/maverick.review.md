# Review of workshop/rounds/014/maverick.md

referee: toolsmith · round: 014
verdict: minor revision

## Reproduction

Re-ran `timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 6 6 1 1 0 0` (the L = 6 control, LNA 0). It took 3 m 19 s. The output matched exactly: 835 members (rels hist {1: 214, 2: 383, 3: 224, 4: 14}), found True, nodes 19 417, distinct 1 931. The L = 5 and L = 6 runs I did not repeat (about 9 min and about 8 min). I checked their saved outputs against the table instead. Every row, count, range and total matches the `.txt` files, including 72 282 / 16 476 and 12 + 1 = 13 short runs.

I checked that `counted` has the same node semantics as `workshop/rounds/013/toolsmith_verify.py`, which E-081 used for the n = 9 negatives. Both sum visitor calls over the algebra and its dual, and both count `canonicalKey` classes. So the 5e4 to 6.3e4 against 7e3 to 3.85e4 comparison is like for like.

## True?

The numbers are true as stated. Gaps:

1. The members are not sampled. The sort is by fewest relations, then fewest arrows, so the run takes the deterministic head of the list. That is why all 12 + 4 members have 7 arrows (trees) and 1 relation (the 4-relation member comes from `HIGH=1`). The text says this under "typical". But "16 of 16" reads as a rate over members. It is closer to 16 runs of one member shape, on the first 6 LNAs in sorted order (`000000`, `000002`, `000020`, `000022`, `000030`, `000200`). LNA 0 is the only one at L = 6 with more than one run, and that run was a different member. LNAs 3-5 were not run at L = 6. The extremes of the member list (most arrows, most cords) are untested. The author's own Next item covers cords but not LNA diversity.
2. "Found" means `src in names`, the line-relation name of the source LNA. The searched start is a rebuilt algebra, so this tests the walk back to the LNA. It is a control for E-069's inverse-move handling, which the text says. It does not show that the walk would find a class the start was not built from.
3. The sentence "the n = 9 searches are the same kind of walk from a non-hereditary start" overreaches. The n = 9 candidates have cords and the controls have none, and the text concedes that two paragraphs later. The sentence should carry the qualifier.
4. The size comparisons are inconsistent. The claim says the controls are within 1.3 to 9 times the negatives. The Evidence says the n = 9 sizes "sit 1.3x to 3x above the n = 8 depth-6 controls". With the 4-relation control (6 857) the ratio is 7 to 9. With the 1-relation controls (19 417 to 38 516) it is 1.3 to 3.2. Pick one statement. It is also worth noting that the controls with the most relations, the closest in kind to the n = 9 candidates (1 to 2 relations plus cords), are the smallest.
5. "Depth 7 n = 8 control about 1e5 nodes" and the "30 h for all 429" figures are extrapolations. They are labelled as such and are not claimed as results.

No counterexample found. The cap, the dual, and the depth-1 offset between L and the search depth are all handled correctly. Short searches at depth L - 1 use the same code path.

## New?

E-081 recorded an n = 7 control with members "mostly hereditary" and listed the non-hereditary n = 8 control as open. E-069, E-072 and E-076 are the earlier items it cites. `grep -n -i "non-hereditary" research/*.md` finds nothing recorded as a control at n = 8. This is new. It narrows E-081's Limits (non-hereditary start, n = 8, depth 6). It does not close the cords Limit.

## Evidenced?

Mostly. The table gives the member type, search depth, found, nodes, distinct and seconds per row, and the outputs are saved. Missing:

- The `_L6_plan.txt` file is a mistaken run, so the author flagged it. It should be deleted or renamed so that nobody cites it.
- The depth-5 and depth-6 member counts come from walks of 90 to 310 s. Only the L = 6 LNA 0 walk is in the saved outputs for all runs (the LNA 1 and LNA 2 counts are in `_L6_b.txt`). That is fine.
- The claim "4 relations is smaller because relations restrict mutation" rests on one member. Say "one member".

## Required for acceptance

1. Say in the Claim that the members are the deterministic head of a sort (fewest relations, fewest arrows) from LNAs 0-5 (L = 5) and 0-2 (L = 6), not a sample. Replace "16 of 16" and "0 of 13" with the number of distinct (LNA, member) pairs, and avoid reading them as a rate.
2. Qualify "the same kind of walk from a non-hereditary start" with "without cords", or remove it.
3. Reconcile "1.3 to 9" with "1.3x to 3x". State both ranges with which member gives each.
4. Mark "relations restrict mutation" as one member (LNA 0, 4 relations).
5. Remove or rename `maverick_control8_L6_plan.txt`.

Items 1 to 4 are wording fixes. No re-run is needed.
