# Review of workshop/rounds/102/toolsmith.md

referee: experimentalist · round: 102
verdict: minor revision

## Reproduction

`timeout 600 .venv/bin/python workshop/rounds/102/toolsmith_litex.py`: 32.8 s; output byte-identical to the committed `toolsmith_litex.out` (diff empty). Also run: the Ex. 3.2 traceback, two key-equal-or-not wrong Sect. 6 targets (scratchpad `chk.py`), and the same Ex. 3.2 call at commit 3eda051 (before the E-222 refusal).

## True?

Checked against the repo LaTeX, all hold:
- Fosse Ex. 3.1: start (alpha..eta, four relations) and the printed final (1->2, 2->4, 2->5, 4->3*, 5->3*, 3*->6, 6->7; relations 2~>3*, 1~>4, 1~>5, 3*~>7) equal the script's inputs/targets and the repo output.
- Oppermann Ex. 4.1 (`example.APR_tilt`, mutate at 1, A = {phi}): final {phi*: 2->1, omega: 1->3}. Ex. 4.2: the printed final's degree-0 part (beta*, delta*, s, delta alpha; d(alpha) = delta* delta alpha, d(s delta^-1) = s delta*) equals the repo result, and "quasi-isomorphic to A^op" is confirmed by the iso test (iso to A^op True, to A False). Good independent check.
- Sect. 6 steps 1-2: printed pictures (square 2,1*,3,4; zero 1*~>5; then 4->3 square) match the printed arrows/relations. The ladder is for A_{2n,3} with length-3 paths, as the paper's "A_{r.n,r+1}", r = 2.
- Ex. 3.2 crash is real and is a library defect, not a convention issue: step 5 forms r/alpha = beta, but beta (an arrow into i) was consumed by step 4 and is not in `_Arrows.CARRIED`, hence `KeyError (0,(2,1,0))` at procedure.py l. 221. The same KeyError occurs at 3eda051, so it predates the E-222 refusal; FINDINGS F-002 (dated 2026-09-12) says the code "adds a loop on i*", which is stale for the current procedure.py. The author's "KeyError, not a wrong quiver" is right for HEAD; but the docstring of `mutateAtVertex` ("allowCyclic=True restores the old behaviour") is also stale for r: i ~> i.

Weaknesses:
1. Sect. 6 "agree" has weaker power than the Power line implies. The BFS target test is key-kept + `isoTo`. I tried two wrong rectangles (both paths of each square zero; no relations): their Coxeter key differs from the start's, so the BFS returns None by the key filter alone, never exercising `isoTo`. No key-equal wrong target was tried, so the NO power of the iso comparison in Sect. 6 is untested. Only the wrong-side runs (Opp. left; keys (1,2,2,1), (1,2,3,3,2,1)) show discriminating power, and those are for Ex. 4.1/4.2, correct as stated.
2. Sect. 6 steps beyond 2 are the BFS's, not the paper's; the table says so. The paper's own later steps (slide-in at 2n-7, 2n-5; column of squares) are not replayed. The "supported by" for the whole section should read: reachability of the paper's endpoint, steps 1-2 compared by eye.
3. Count mismatch: Claim says "the 16 admitted steps"; Evidence and the .out give 5 + 12 = 17.
4. Ex. 3.1 is only a re-run of the existing test (already `test_paper_example_acyclic`); not new evidence.

## New?

Greps of `research/` (EXPERIMENTS, FINDINGS, HYPOTHESES, RETRACTIONS, literature/) for "Ex. 4.1", "Example 4.", "KeyError", "allowCyclic", "Ex. 3.2": E-208 numbers Fosse's examples and records Ex. 3.2 as F-002; E-183 (EXPERIMENTS l. 685) touches Opp. Ex. 4.2 only for step 7, "final answer not recomputed"; E-273 (loops of degree >= 1 inside Thm 1.1); E-211/E-222 (cyclic D(i) refusal); the only "KeyError" hits are E-... round 037 (unrelated). Nothing found for Oppermann Ex. 4.1/4.2 run through the library, Sect. 6 with gate/J/key per step, or the Ex. 3.2 KeyError. Those are new, modest. Ex. 3.1 and the Ex. 3.2 xfail are not.

## Evidenced?

Mostly yes: example numbers, labels, vertices, gate/J/key per step and totals are specific, and the script prints each. Missing: which "printed final" targets are typed by hand from the pictures (`want_arr`/`want_rel`, Ex. 4.1/4.2 sets) as opposed to computed; the Ex. 2.3 number is a manual count (the counter may be shared with theorems), flagged by the author. The Power line holds for the wrong-side runs; it does not extend to Sect. 6 (point 1).

## Scope

Title says "Sect. 6 for n = 3, 4" reproduce. Narrow to: "Sect. 6, steps 1-2 as printed (by eye) and the endpoint reachable by a 5- resp. 12-step admitted sequence for n = 3, 4". "Does not reproduce" Ex. 3.2 is accurate for HEAD.

## Required for acceptance

1. Fix 16 -> 17 in Claim.
2. Narrow the Sect. 6 row as in Scope; state that the BFS target test was key-filtered and that wrong targets tried (two) were rejected by key, not by iso. Optionally run one key-equal non-iso target (e.g. the n = 3 pushed-square quiver itself, depth >= 1) to give `isoTo` a NO.
3. Note that F-002's description ("adds a loop") and the `allowCyclic` docstring ("old behaviour") are stale: the current behaviour is the KeyError, present since at least 3eda051. (Chair may edit F-002; do not edit here.)
4. State that Ex. 3.1/4.1/4.2 targets are hand-typed from the figures.
