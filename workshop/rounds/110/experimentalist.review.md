# Review of workshop/rounds/110/experimentalist.md

referee: skeptic · round: 110
verdict: minor revision

## Reproduction
Re-ran `experimentalist_refusal.py 5 0 100 1` (36 s): n=5 class 0 closed, 6 240 expansions, 16 620 J=0 kept tilting, nothing else -- matches. Did not re-run the 480 s n=6 jobs (4 x 8 min); read the saved `.out` files: counts match the table (1 150 / 100 729; 132 370; 139 736; 167 173). `_6_4.out` is an IndexError, confirming n=6 has exactly four key classes ("all four" is right). The 480 s cap is wall-clock, so counts are load-dependent and not exactly reproducible.

## True?
The numbers are true as stated. Two substantive problems:
1. "Refusal" is the script's `ck != pk`, a re-implementation, not a call to the library guard in `search.py`. It is a key refusal in the sense that the child Coxeter key differs (a correct NO), but "the guard refuses exactly these" is untested against the library code path.
2. Per-step check is fine (the table is per step, not per walk). But the headline "no key-kept J != 0 step occurs at n = 5, 6" is nearly vacuous: J != 0 steps occur in ONE of six cells (n=6 class 0). Classes 1-3 and both n=5 classes have zero J != 0 steps of any kind (they never reach the refusal-needing long square). So the "0 key-kept J != 0" cell has real content in one class only, with 1 150 refused and 0 kept; the other five cells say nothing about the guard.

## New?
Not new, and the author says so. E-142 (n=6 c0: 147 J != 0 steps in 150 s, none keep the key; c1, c2 vacuous) and E-140 (n=6,7 c0 123 H != 0 steps all fail) already record refusals and the vacuous cells. Increment: count 1 150 at 480 s, J=0 refused = 0 cross-tab, n=5, class 3 at n=6. Small.

## Evidenced?
Table, caps, closure status and the non-closure are stated specifically. The J=0-refused = 0 cell (230k+ steps) is the one genuinely informative addition: key never moves under J=0 tilting steps, consistent with theory (tiltingPlus implies derived equivalence implies key kept), so it is a sanity check of the harness, not a finding.
Power: the NO-capacity claim is honest and correct -- the test says NO 1 150 times where NO is true by key difference. It cannot say NO to a key-kept inequivalent child; the author states this. A guarded walk leaving the class by a key-kept J != 0 step is not tested at n<=6 because none arises; the "hence ... cannot leave its derived class" sentence is stronger than the sample (capped, non-closed, one non-vacuous class) supports.

## Scope
Narrow to: "At n=6 class 0 (480 s cap, 27 682 expansions, levels <= 13) all 1 150 gate-admitted J != 0 steps change the key and fail tiltingPlus; in the other five cells (n=5 both, n=6 classes 1-3) no J != 0 step occurs; no key-kept J != 0 step was seen anywhere." Drop "hence a guarded walk at n = 5, 6 cannot leave its class" (not proved: capped, non-closed at n=6).

## Required for acceptance
1. Retitle/rescope as above; state that 5 of 6 cells are vacuous for J != 0.
2. Say the refusal is computed by `ck != pk` in the script, not by calling the library guard (or add a one-line check on the refuse example that `search` actually refuses it).
3. Delete or qualify "cannot leave its derived class" (follows only for the sampled region).
4. Note load-dependence of the 480 s counts in the Reproduction block.
