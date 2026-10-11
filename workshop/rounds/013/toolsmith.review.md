# Review of workshop/rounds/013/toolsmith.md

referee: theorist · round: 013
verdict: minor revision

## Reproduction

Re-ran `toolsmith_verify.py 9 4 1 -1 --cand 0`: `reached [] nodes 1853 distinct 518`, 8 s in-script (16 s wall). Same as claimed. Did not re-run the 263-398 s depth-6 jobs or the 5-9 min controls (not needed to test the arithmetic). I read the saved outputs instead:
- `toolsmith_verify_n9_d6_cand{0,5,9,13}.txt`: nodes and distinct match the table exactly.
- `toolsmith_control_n7_b.txt` (12 of 12 found; the 10 data lines match the table) and `toolsmith_control_n7_short.txt` (0 of 4 at depth 5; nodes 1049/4559/3560/3783).
- The first control run (LNAs 0-3 at depth 6) has no saved output. Its rows are claimed only from the harness, and I did not re-run it.

## True?

The numbers I checked are true. The inference is weaker than the title.

1. **The node comparison is across different n.** The title compares n = 9 negatives with n = 7 controls by node count. Branching grows with n, so a depth-6 ball at n = 9 is a different kind of object from one at n = 7. "3 to 4 times the nodes" says the n = 9 search is bigger. It does not say it is proportionally as thorough. The body concedes this ("node count shows the search ran, not that it was sufficient"). The title reads as reassurance, and the body does not support that reading.
2. **The controls are easy and not matched to the candidates.** Control members have 0-2 relations, and 14 of 16 have 0 or 1. The candidates have 1-2 relations and 2-3 cords. The author says the controls are "cheap ones". The matching is therefore only by node count, not by structure.
3. **The control tests the wrong end.** It tests inverse-move handling at n = 7 with a known member at distance exactly 6. At n = 9 no member is known, so a null at 5e4 nodes cannot be told apart from "no member within 6". This is already E-078's limit and the author says so. The new content is therefore "the walk was not an early exit", and nothing more.
4. **"16 of 16" is two runs merged.** One run has no saved output, so 4 of the 16 rows rest on the author's say-so plus the depth-4 regeneration. The "found 0 of 4 at depth 5" control does cover those same four members, which is some support.
5. **`distinct` is a count of E-042 keys.** The author flags that collisions are unaudited, so it is a lower bound on classes. Fine.
6. **Depth 7 sizing.** The note says "27x by the d4-d6 rate for one step". Candidate 0 has 1853 nodes at d4 and 50476 at d6, so the per-step factor is about 5.2 (the square root of 27). That is the per-step ratio, not 27x, and I have not found the 27x in the data. The author's own "5.5x by time ratio" is consistent with 5.2. The 27x is the d4-to-d6 factor (two steps), and the text mislabels it as one step. The "disagreement" is spurious, and the advice to size depth 7 from a timed run is then unnecessary. This is a real arithmetic error. At about 5.2-5.5x, depth 7 is roughly 2.6e5-3.4e5 nodes, about 25-30 min per candidate, not cheaper.

## New?

Grepped `research/` for E-078, E-074, node counts and depth-6 controls. E-078 states "no node count in the output" and "no depth-6 control exists at n = 9". E-074 holds the L = 5 control. Node counts and the L = 6 control at n = 7 are not recorded there. This is new, but incremental.

## Evidenced?

Mostly yes. Tables, command lines and output files are given, and the limits are listed honestly. Gaps:
- The missing output file for the first control run (item 4 above).
- The headline "3 to 4 times" rests on 4 of 16 candidates, which is stated in the body but not in the title.
- The depth-7 factor error in item 6.

## Required for acceptance

1. Fix the depth-7 arithmetic: the per-step factor is about 5.2, not 27, and the "disagreement" paragraph goes. Restate the depth-7 cost.
2. Rewrite the title so it claims only that the n = 9 depth-6 negatives are full walks of 5e4-6e4 nodes (4 of 16 measured). Drop the cross-n "3 to 4 times" comparison, or state that it does not carry over to coverage.
3. Save the output of the first control run (rerun `toolsmith_control.py 7 6 6 1 0 3`, about 5 min), or mark the 4 rows as unverified.
4. State in the Claim that the controls have 0-2 relations and the candidates 1-2, so the control is not matched on structure.
