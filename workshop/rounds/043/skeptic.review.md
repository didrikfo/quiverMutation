# Review of workshop/rounds/043/skeptic.md

referee: theorist · round: 043
verdict: minor revision

## Reproduction

- `skeptic_c2.py 7 1 2 500 30 guard` and `... 7 2 3 500 30 guard` (500 s each, run in parallel). Output matched the table to the count. c1: 75 steps, 67 distinct, D=0 in 12 H1/H2 plus 1 non-H1/H2 (=13), H1/H2 low terms (1,-1) 2, (2,1) 48, (3,2) 4, orbit s=1 in 22 and none in 44. c2: 64 distinct, D=0 in 9, (3,1) in 52, s=10 in 9, none in 52. The 12 and 9 are strict H1/H2.
- `skeptic_zero.py skeptic_zero_n7c2.pkl` (about 30 s). The 4 pickled c2 steps rebuild from arrows and relations. For each: gate True, tiltingPlus False, key(parent) = key(child) = hand-built child key = (1,1,0,-2,-2,0,1,1), canonical key of the library child equals the hand-built child, Cartan incongruent, dim 22/17 or 19/15.
- Not re-run: n=6 guard-off (540 s), n=7 c0, the random m=5/6/7 runs, and the tguard runs. The tguard c1/c2 counts, the random-parent claims (3), (4) and the c2 D=0 "s=10" cell are therefore unchecked here.

## True?

Claim (2) holds. The three things I was asked to check:

- **Child key computed on the child.** Yes. `search._coxeterKeyOrNone(ch)` is applied to the reduced child `ch`, and the hand rebuild confirms it independently.
- **Gate-admitted.** Yes. `mutationIsPossibleAtVertex` is checked before the step is counted, and the rebuild prints gate True.
- **Key of E-142.** The same function is used as in the class walk. The D=0 test is equivalent to key equality here, because det C = 1 and det(xC+C^T) is then the Coxeter polynomial. One residual caveat: "class 1/2" are indices into classes sorted by (size, str(key)), with sizes [8,12,14,58,64,108] at n=7. I did not check that this ordering matches the class numbering of E-142 (`experimentalist_keyoff.py`). The report should print the key per class (c2 is the one shown, (1,1,0,-2,-2,0,1,1)).

Gaps, none fatal:

1. **Not an LNA walk.** These are not steps from LNA parents. The parents are walk descendants, and every step has tiltingPlus False. In the tguard runs the walk only follows steps with tiltingPlus True, so those parents are tilting-reachable from LNAs, but the D=0 step itself is never tilting. This does contradict E-142's reading, because E-142 also walked descendants (a guarded walk, not just LNA parents). It does not contradict E-142's data: E-142 had 2 J != 0 steps in c1 and none in c2, at 150 s and shallower depth. The report's wording "depth/time artefact" is an inference. It was not tested by re-running E-142's exact command at 500 s.
2. **Not derived equivalence.** Key-equal does not mean derived equivalent. The report says this (Cartan incongruent). It should say outright that the key guard therefore does not separate derived classes off J=0 steps, so it is no longer evidence for H-015 outside class 0.
3. **c1 D=0 steps are not independently rebuilt.** The pickle fails to rebuild (parallel arrows). The c1 count rests on one script. The c2 steps are independently confirmed.
4. **Claims (3) and (4) go beyond the evidence.** The "orbit relation fails" and "c_2 != 0 realisable" statements come from non-LNA parents (random acyclic m=5..7). The report says so, but the title and the (4) "none has D=0" line imply more than the random runs show. Claim (3) for c1/c2 is fine: s is searched only over |s| <= 40 with float `allclose`, so "no s" really means "no s <= 40 numerically". For F of finite order this is exhaustive, but it should be stated.

## New?

Grepped FINDINGS, HYPOTHESES and RETRACTIONS for "key-preserving" and "orbit relation": nothing. EXPERIMENTS.md has E-142 (none key-preserving, c1 2 steps, c2 vacuous). That is contradicted in reading, not in data. E-145 line (EXPERIMENTS.md, the P1/P2 entry) covers H1/H2, the orbit relation s=1 or -2, and class 0, and its own text says the s=-2 case has no argument. It also says "50 off-shape steps (u = e_a+e_b-e_i, or H1 fails) also have lowest term x^2". E-142's control found 55 key-preserving J != 0 steps among random acyclic parents, so key-preserving J != 0 steps were known off the LNA class. New: they occur inside the LNA key classes 1 and 2, with H1/H2, and the orbit relation fails there.

## Evidenced?

Mostly. The counts are specific, and the c2 example is given in full and rebuilds. Missing:

- the key printed for each class, and the class ordering used;
- the E-142 command re-run at matched budget, which would show the change in count is due to time/depth;
- the "21 steps where library key and D agree" table, which is cited but not shown;
- the full list of D=0 hits. Only the first 4 per run are pickled.

## Required for acceptance

1. Print the class key for n=7 index 1 and 2, and confirm that index matches E-142's c1 and c2.
2. Run `experimentalist_keyoff.py 7 1 ...` (E-142's command) at 500 s and report whether it now finds the D=0 steps, to support the "depth/time artefact" statement.
3. Fix the c1 rebuild (parallel arrows), or state that the c1 count is single-script.
4. Narrow the title: say "off LNA parents, in walk descendants, tiltingPlus False at the step", not just "key-preserving steps exist". State the range |s| <= 40.
5. Mark E-142's c1 and c2 cells in the Next section as non-vacuous (done), and add that the E-140/E-143 "J != 0 moves the key" law is false at n=7 c1, c2 on walk descendants (not only on random parents).
