# Review of workshop/rounds/047/toolsmith.md

referee: theorist · round: 047
verdict: minor revision

## Reproduction

`toolsmith_witness.py --all`: 19 s, output matches the table. 10 pairs, each halves 3+3, total 6, one tied meeting each. Paths are [1,1,2]/[1,2,1] for the member-start pairs and [-8,-8,-7]/[-8,-7,-8] for the dual-start pairs (the submission writes the dual ones as 8,8,7). TOTAL edges 60, gate 60, J0 60, key 60, and the replay ends on the meeting key in 10 of 10. `pytest tests/test_merges_witness.py tests/test_merge_decisions.py -m "not slow"`: 27 passed. I did not run the `merges.py 10 --depths 5 --witness` line (the submission says it is untested by data).

## True?

Nothing wrong in what is stated.
- The "no witness of length <= 5" step holds: depth <= 3 per side covers every split of a total <= 5, and the dual option is searched on both sides.
- It holds only inside the key-guarded gate graph, as claim (2) says.
- Claim (3) is not what the title suggests. The edges tested are the A-half forward and the B-half forward from B. The route A to B traverses the B-half inverted, so "all 60 edges pass tiltingPlus" does not show that the 6-step A-to-B route is a tilting path. The body states this caveat; the title does not.
- "Shortest" is among the first paths stored per key (one path per key, ties not enumerated). Disclosed.
- The script uses `workshop/rounds/001/scholar_h015.py` for `tiltingPlus` via exec. I did not audit that, but E-148 uses the same J = 0 test.

## New?

Mostly not. E-148 already records: the same 10 pairs, meeting at total 6 with split (3,3), and all 2396 gate-admitted edges J = 0. F-041 records the meeting depth 3+3. E-147 is the source of the reverse-step caveat. HYPOTHESES (H-015 ledger, round 045) lists the deep-merge replay as the thing still to close. Genuinely new: the explicit paths, the per-edge check on the witness itself (60 edges), the opt-in `--witness` storage in `merges.py`, and the statement that the B-half direction is untested. The author says as much in "Prior record".

## Evidenced?

Yes for the script result: the table, the command and the run time are given and they reproduce. Two gaps:
- The `witnesses` filtering in `merges.py main` has no data test. It is disclosed, and the n = 5 unit test covers the search layer only.
- The "one local move" reading (two path shapes) is a pattern in 10 rows, not a derivation. It is correctly labelled as agreeing with E-148 and not as a proof.

## Scope

n = 8, these 10 pairs, depth 3 per side, guarded graph with duals. The claim text matches this. The title should carry the same limits: add "(key-guarded graph; edges as searched, B-half not inverted)" or similar.

## Required for acceptance

1. Narrow the title: "all 60 edges pass tiltingPlus" should say the B-half edges are tested forward from B, not along the A-to-B route.
2. Say in the claim, not only the body, that "no witness of length <= 5" is relative to the key-guarded gate graph.
3. Optional and small: run `replay` on the inverse B-half (or state it is left to the Skeptic, as already written under Next). If done in one sitting, report the J = 0 count over the 30 inverse edges. If not, mark it `[next round]`.
